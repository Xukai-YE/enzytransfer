#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Enhanced mutation retry with automatic sequence fetching.

NEW FEATURES:
1. Auto-fetch wildtype_sequence from UniProt when missing
2. Try canonical sequence first, then isoforms (1-6)
3. Local disk cache to avoid redundant API calls
4. Rate limiting with exponential backoff retry

This script processes YAML files from the issues directory and attempts to:
1. Cluster files by UniProt ID (organized into separate folders)
2. Auto-fetch missing wildtype sequences from UniProt
3. Build regex patterns from mutation positions to locate the correct region in the sequence
4. Calculate position offset and re-apply mutations with corrected positions

Features:
- Checkpoint support: resume from last interruption
- Detailed progress output
- Organized output by UniProt ID
- **NEW**: Automatic sequence retrieval with caching

Usage:
    python step1_enhanced.py --input ./mutation_exports/issues
    
    # Resume from checkpoint
    python step1_enhanced.py --input ./mutation_exports/issues --resume
    
    # Custom cache directory
    python step1_enhanced.py --input ./mutation_exports/issues --cache-dir ./uniprot_cache
"""

import argparse
import csv
import glob
import json
import os
import re
import sys
import time
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any

import yaml


# ==================== UniProt Sequence Cache System ====================

class UniProtCache:
    """
    Disk-based cache for UniProt sequences with rate limiting.
    
    Features:
    - Persistent cache on disk (JSON format)
    - Rate limiting (1 request per second by default)
    - Exponential backoff retry
    - Separate caching for canonical and isoform sequences
    """
    
    def __init__(self, cache_dir: str = "./uniprot_cache", rate_limit: float = 1.0):
        """
        Args:
            cache_dir: Directory to store cache files
            rate_limit: Minimum seconds between API requests
        """
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.cache_file = self.cache_dir / "sequence_cache.json"
        self.rate_limit = rate_limit
        self.last_request_time = 0
        
        # Load existing cache
        self.cache = self._load_cache()
        
        # Statistics
        self.stats = {
            'cache_hits': 0,
            'cache_misses': 0,
            'api_requests': 0,
            'api_errors': 0
        }
    
    def _load_cache(self) -> Dict:
        """Load cache from disk."""
        if self.cache_file.exists():
            try:
                with open(self.cache_file, 'r', encoding='utf-8') as f:
                    return json.load(f)
            except Exception as e:
                print(f"[WARNING] Failed to load cache: {e}")
        return {}
    
    def _save_cache(self):
        """Save cache to disk."""
        try:
            with open(self.cache_file, 'w', encoding='utf-8') as f:
                json.dump(self.cache, f, indent=2, ensure_ascii=False)
        except Exception as e:
            print(f"[WARNING] Failed to save cache: {e}")
    
    def _wait_for_rate_limit(self):
        """Enforce rate limiting between API requests."""
        elapsed = time.time() - self.last_request_time
        if elapsed < self.rate_limit:
            sleep_time = self.rate_limit - elapsed
            time.sleep(sleep_time)
        self.last_request_time = time.time()
    
    def _make_cache_key(self, accession: str, isoform: Optional[int] = None) -> str:
        """Generate cache key for a sequence."""
        if isoform is None:
            return f"{accession}:canonical"
        return f"{accession}:isoform-{isoform}"
    
    def get_cached_sequence(self, accession: str, isoform: Optional[int] = None) -> Optional[str]:
        """
        Get sequence from cache if available.
        
        Args:
            accession: UniProt accession
            isoform: Isoform number (None for canonical)
        
        Returns:
            Cached sequence or None if not in cache
        """
        cache_key = self._make_cache_key(accession, isoform)
        
        if cache_key in self.cache:
            self.stats['cache_hits'] += 1
            return self.cache[cache_key]
        
        self.stats['cache_misses'] += 1
        return None
    
    def store_sequence(self, accession: str, sequence: str, isoform: Optional[int] = None):
        """
        Store sequence in cache.
        
        Args:
            accession: UniProt accession
            sequence: Amino acid sequence
            isoform: Isoform number (None for canonical)
        """
        cache_key = self._make_cache_key(accession, isoform)
        self.cache[cache_key] = sequence
        self._save_cache()
    
    def fetch_canonical_sequence(self, accession: str, max_retries: int = 3) -> Optional[str]:
        """
        Fetch canonical sequence from UniProt with retry logic.
        
        Args:
            accession: UniProt accession
            max_retries: Maximum number of retry attempts
        
        Returns:
            Sequence string or None if fetch failed
        """
        import requests
        
        # Check cache first
        cached = self.get_cached_sequence(accession, isoform=None)
        if cached:
            return cached
        
        # Fetch from API with exponential backoff
        url = f"https://rest.uniprot.org/uniprotkb/{accession}.fasta"
        
        for attempt in range(max_retries):
            try:
                # Rate limiting
                self._wait_for_rate_limit()
                
                # Make request
                self.stats['api_requests'] += 1
                response = requests.get(url, timeout=30)
                response.raise_for_status()
                
                # Parse FASTA
                text = response.text.strip()
                if not text.startswith(">"):
                    raise ValueError(f"Invalid FASTA format for {accession}")
                
                # Extract sequence (skip header line)
                lines = text.splitlines()
                sequence = "".join(line.strip() for line in lines[1:] if line and not line.startswith(">"))
                
                if not sequence:
                    raise ValueError(f"Empty sequence for {accession}")
                
                # Cache and return
                self.store_sequence(accession, sequence, isoform=None)
                return sequence
            
            except requests.exceptions.HTTPError as e:
                if e.response.status_code == 404:
                    # Not found - don't retry
                    self.stats['api_errors'] += 1
                    return None
                
                # Other HTTP errors - retry with backoff
                if attempt < max_retries - 1:
                    wait_time = 2 ** attempt  # Exponential backoff: 1s, 2s, 4s
                    print(f"[WARNING] HTTP error fetching {accession} (attempt {attempt + 1}/{max_retries}): {e}")
                    print(f"          Retrying in {wait_time}s...")
                    time.sleep(wait_time)
                else:
                    self.stats['api_errors'] += 1
                    print(f"[ERROR] Failed to fetch {accession} after {max_retries} attempts: {e}")
                    return None
            
            except Exception as e:
                # Other errors - retry with backoff
                if attempt < max_retries - 1:
                    wait_time = 2 ** attempt
                    print(f"[WARNING] Error fetching {accession} (attempt {attempt + 1}/{max_retries}): {e}")
                    print(f"          Retrying in {wait_time}s...")
                    time.sleep(wait_time)
                else:
                    self.stats['api_errors'] += 1
                    print(f"[ERROR] Failed to fetch {accession} after {max_retries} attempts: {e}")
                    return None
        
        return None
    
    def fetch_isoform_sequence(self, accession: str, isoform: int, max_retries: int = 3) -> Optional[str]:
        """
        Fetch isoform sequence from UniProt with retry logic.
        
        Args:
            accession: UniProt accession
            isoform: Isoform number
            max_retries: Maximum number of retry attempts
        
        Returns:
            Sequence string or None if fetch failed
        """
        import requests
        
        # Check cache first
        cached = self.get_cached_sequence(accession, isoform=isoform)
        if cached:
            return cached
        
        # Fetch from API with exponential backoff
        query = f"accession:{accession}-{isoform}"
        url = f"https://rest.uniprot.org/uniprotkb/stream?format=fasta&query={query}"
        
        for attempt in range(max_retries):
            try:
                # Rate limiting
                self._wait_for_rate_limit()
                
                # Make request
                self.stats['api_requests'] += 1
                response = requests.get(url, timeout=30)
                response.raise_for_status()
                
                # Parse FASTA
                text = response.text.strip()
                if not text or not text.startswith(">"):
                    # Isoform doesn't exist
                    return None
                
                # Extract sequence
                lines = text.splitlines()
                sequence = "".join(line.strip() for line in lines[1:] if line and not line.startswith(">"))
                
                if not sequence:
                    return None
                
                # Cache and return
                self.store_sequence(accession, sequence, isoform=isoform)
                return sequence
            
            except requests.exceptions.HTTPError as e:
                if e.response.status_code == 404:
                    # Not found - don't retry
                    return None
                
                # Other HTTP errors - retry with backoff
                if attempt < max_retries - 1:
                    wait_time = 2 ** attempt
                    print(f"[WARNING] HTTP error fetching {accession}-{isoform} (attempt {attempt + 1}/{max_retries}): {e}")
                    print(f"          Retrying in {wait_time}s...")
                    time.sleep(wait_time)
                else:
                    print(f"[ERROR] Failed to fetch {accession}-{isoform} after {max_retries} attempts: {e}")
                    return None
            
            except Exception as e:
                # Other errors - retry with backoff
                if attempt < max_retries - 1:
                    wait_time = 2 ** attempt
                    print(f"[WARNING] Error fetching {accession}-{isoform} (attempt {attempt + 1}/{max_retries}): {e}")
                    print(f"          Retrying in {wait_time}s...")
                    time.sleep(wait_time)
                else:
                    print(f"[ERROR] Failed to fetch {accession}-{isoform} after {max_retries} attempts: {e}")
                    return None
        
        return None
    
    def fetch_all_sequences(self, accession: str, max_isoforms: int = 6) -> Dict[str, str]:
        """
        Fetch canonical and all isoforms for a UniProt accession.
        
        Args:
            accession: UniProt accession
            max_isoforms: Maximum number of isoforms to try
        
        Returns:
            Dict mapping 'canonical' or 'isoform-N' to sequence
        """
        sequences = {}
        
        # Fetch canonical
        canonical = self.fetch_canonical_sequence(accession)
        if canonical:
            sequences['canonical'] = canonical
        
        # Fetch isoforms
        for iso_num in range(1, max_isoforms + 1):
            isoform = self.fetch_isoform_sequence(accession, iso_num)
            if isoform:
                sequences[f'isoform-{iso_num}'] = isoform
        
        return sequences
    
    def print_stats(self):
        """Print cache statistics."""
        total_queries = self.stats['cache_hits'] + self.stats['cache_misses']
        hit_rate = 100 * self.stats['cache_hits'] / total_queries if total_queries > 0 else 0
        
        print("\n📊 UniProt Cache Statistics:")
        print(f"  Cache hits:    {self.stats['cache_hits']:4d}")
        print(f"  Cache misses:  {self.stats['cache_misses']:4d}")
        print(f"  Hit rate:      {hit_rate:5.1f}%")
        print(f"  API requests:  {self.stats['api_requests']:4d}")
        print(f"  API errors:    {self.stats['api_errors']:4d}")
        print(f"  Cache file:    {self.cache_file}")


# ==================== Original Code (with enhancements) ====================

# -------- Checkpoint manager --------
class CheckpointManager:
    """Manages checkpoint data for resumable processing."""
    
    def __init__(self, checkpoint_file: str):
        self.checkpoint_file = checkpoint_file
        self.data = self._load()
    
    def _load(self) -> Dict:
        if os.path.exists(self.checkpoint_file):
            try:
                with open(self.checkpoint_file, 'r') as f:
                    return json.load(f)
            except Exception:
                pass
        return {
            'processed_files': [],
            'failed_files': [],
            'stats': {
                'total': 0,
                'success': 0,
                'failed': 0,
                'skipped': 0,
                'auto_fetched': 0  # NEW: Track auto-fetched sequences
            },
            'last_update': None
        }
    
    def save(self):
        self.data['last_update'] = datetime.now().isoformat()
        with open(self.checkpoint_file, 'w') as f:
            json.dump(self.data, f, indent=2)
    
    def is_processed(self, file_path: str) -> bool:
        return file_path in self.data['processed_files']
    
    def mark_processed(self, file_path: str, success: bool = True):
        if file_path not in self.data['processed_files']:
            self.data['processed_files'].append(file_path)
        if not success and file_path not in self.data['failed_files']:
            self.data['failed_files'].append(file_path)
    
    def update_stats(self, **kwargs):
        for key, value in kwargs.items():
            if key in self.data['stats']:
                self.data['stats'][key] += value
    
    def get_stats(self) -> Dict:
        return self.data['stats']
    
    def clear(self):
        """Clear checkpoint data."""
        self.data = {
            'processed_files': [],
            'failed_files': [],
            'stats': {'total': 0, 'success': 0, 'failed': 0, 'skipped': 0, 'auto_fetched': 0},
            'last_update': None
        }
        if os.path.exists(self.checkpoint_file):
            os.remove(self.checkpoint_file)


# -------- Progress display helpers --------
def print_header(text: str):
    """Print a formatted header."""
    line = "=" * 70
    print(f"\n{line}")
    print(f"  {text}")
    print(f"{line}\n")


def print_progress(current: int, total: int, prefix: str = "", suffix: str = ""):
    """Print progress bar."""
    bar_length = 40
    filled_length = int(bar_length * current / total) if total > 0 else 0
    bar = '█' * filled_length + '-' * (bar_length - filled_length)
    percent = f"{100 * current / total:.1f}%" if total > 0 else "0.0%"
    print(f'\r{prefix} |{bar}| {percent} {suffix}', end='', flush=True)
    if current == total and total > 0:
        print()  # New line on completion


def format_time(seconds: float) -> str:
    """Format seconds into human-readable time."""
    if seconds < 60:
        return f"{seconds:.1f}s"
    elif seconds < 3600:
        return f"{int(seconds // 60)}m {int(seconds % 60)}s"
    else:
        hours = int(seconds // 3600)
        minutes = int((seconds % 3600) // 60)
        return f"{hours}h {minutes}m"


# -------- UniProt ID handling --------
def has_multiple_uniprot_ids(uniprotid: str) -> bool:
    """Check if a uniprotid string contains multiple IDs."""
    if not uniprotid:
        return False
    
    uniprotid = str(uniprotid).strip()
    
    # Check for common separators
    for sep in [';', ',', '|']:
        if sep in uniprotid:
            return True
    
    return False


def get_primary_uniprot_id(uniprotid: str) -> str:
    """Extract the primary (first) UniProt ID from a string that may contain multiple IDs."""
    if not uniprotid:
        return ""
    
    uniprotid = str(uniprotid).strip()
    
    # Split by common separators and take the first ID
    for sep in [';', ',', '|']:
        if sep in uniprotid:
            return uniprotid.split(sep)[0].strip()
    
    return uniprotid


def get_original_uniprot_id(prot: Dict) -> str:
    """Get the original uniprotid from protein dict (may contain multiple IDs)."""
    return str(prot.get("uniprotid") or prot.get("uniprot_id") or prot.get("uniprot") or "")


# -------- Parse variant_description --------
def _tokenize_variants(desc: str) -> List[str]:
    if not desc:
        return []
    s = desc.strip()
    s = re.sub(r'^\s*p\.\[?', '', s, flags=re.IGNORECASE)
    s = re.sub(r'\]$', '', s)
    parts = re.split(r'[;,\s]+', s)
    return [p for p in parts if p]


THREE_TO_ONE = {
    "Ala": "A", "Cys": "C", "Asp": "D", "Glu": "E", "Phe": "F", "Gly": "G", "His": "H",
    "Ile": "I", "Lys": "K", "Leu": "L", "Met": "M", "Asn": "N", "Pro": "P", "Gln": "Q",
    "Arg": "R", "Ser": "S", "Thr": "T", "Val": "V", "Trp": "W", "Tyr": "Y", "Ter": "*",
    "Sec": "U", "Pyl": "O", "Xaa": "X"
}


def three_to_one(aa3: str) -> str:
    key = aa3[:1].upper() + aa3[1:].lower()
    if key not in THREE_TO_ONE:
        raise ValueError(f"Unknown 3-letter AA code: {aa3}")
    return THREE_TO_ONE[key]


def parse_variant_description(desc: Optional[str]) -> List[Tuple[str, int, str]]:
    """
    Returns a list of (fromAA, pos, toAA).
    Note: Single-letter amino acids are automatically converted to uppercase.
    """
    if not desc:
        return []
    tokens = _tokenize_variants(desc)
    muts: List[Tuple[str, int, str]] = []
    for t in tokens:
        if t == "=":
            continue
        
        # Match single-letter format (e.g., F236S or f236s)
        m1 = re.fullmatch(r'([A-Za-z\*])(\d+)([A-Za-z\*])', t)
        if m1:
            # Convert to uppercase for consistency
            from_aa = m1.group(1).upper()
            to_aa = m1.group(3).upper()
            muts.append((from_aa, int(m1.group(2)), to_aa))
            continue
        
        # Match three-letter format (e.g., Phe236Ser)
        m2 = re.fullmatch(r'([A-Za-z]{3})(\d+)([A-Za-z]{3}|=)', t)
        if m2:
            if m2.group(3) == "=":
                continue
            aa0 = three_to_one(m2.group(1))
            aa1 = three_to_one(m2.group(3))
            muts.append((aa0, int(m2.group(2)), aa1))
            continue
        
        raise ValueError(f"Unrecognized mutation token: '{t}' (from '{desc}')")
    
    # de-duplicate
    seen = set()
    out = []
    for a, p, b in muts:
        key = (a, p, b)
        if key not in seen:
            seen.add(key)
            out.append(key)
    return out


# -------- Pattern-based mutation retry --------
def try_single_mutation_with_offset(seq: str, orig_aa: str, pos: int, new_aa: str, 
                                     max_offset: int = 10) -> Tuple[bool, int, str]:
    """Try to apply a single mutation with position offset search."""
    # Try exact position first
    if 1 <= pos <= len(seq) and seq[pos - 1] == orig_aa:
        return (True, pos, f"Found at exact position {pos}")
    
    # Search within offset range
    for offset in range(-max_offset, max_offset + 1):
        if offset == 0:
            continue  # Already tried
        
        actual_pos = pos + offset
        if 1 <= actual_pos <= len(seq) and seq[actual_pos - 1] == orig_aa:
            return (True, actual_pos, f"Found at position {actual_pos} (offset: {offset:+d})")
    
    return (False, 0, f"Amino acid '{orig_aa}' not found near position {pos} (±{max_offset})")


def build_regex_pattern(mutations: List[Tuple[str, int, str]]) -> Tuple[str, int, int]:
    """Build a regex pattern from mutations to locate the region in the sequence."""
    if not mutations:
        raise ValueError("No mutations provided")
    
    # Sort by position
    sorted_muts = sorted(mutations, key=lambda x: x[1])
    min_pos = sorted_muts[0][1]
    max_pos = sorted_muts[-1][1]
    
    # Build pattern: use original amino acids at mutation positions, '.' for gaps
    pattern_parts = []
    current_pos = min_pos
    
    for orig_aa, pos, _ in sorted_muts:
        gap = pos - current_pos
        if gap > 0:
            pattern_parts.append('.' * gap)
        pattern_parts.append(orig_aa)
        current_pos = pos + 1
    
    pattern = ''.join(pattern_parts)
    return pattern, min_pos, max_pos


def find_pattern_in_sequence(seq: str, pattern: str, expected_pos: int, 
                             tolerance: int = 50) -> Optional[int]:
    """Find the pattern in the sequence, searching around the expected position."""
    # Convert to 0-based for searching
    search_start = max(0, expected_pos - 1 - tolerance)
    search_end = min(len(seq), expected_pos - 1 + tolerance + len(pattern))
    search_region = seq[search_start:search_end]
    
    try:
        match = re.search(pattern, search_region)
        if match:
            # Convert back to 1-based position
            return search_start + match.start() + 1
        
        # If not found in tolerance region, search entire sequence
        match = re.search(pattern, seq)
        if match:
            return match.start() + 1
    except re.error:
        pass
    
    return None


def apply_mutations_with_offset(seq: str, mutations: List[Tuple[str, int, str]], 
                                offset: int) -> Tuple[str, List[str]]:
    """Apply mutations with position offset."""
    s = list(seq)
    messages = []
    
    for orig_aa, orig_pos, new_aa in mutations:
        corrected_pos = orig_pos + offset
        
        if not (1 <= corrected_pos <= len(s)):
            msg = f"Position {orig_pos} (corrected to {corrected_pos}) out of range"
            messages.append(msg)
            raise ValueError(msg)
        
        actual_aa = s[corrected_pos - 1]
        
        if actual_aa != orig_aa:
            msg = f"Mismatch at {orig_pos} (corrected to {corrected_pos}): expected '{orig_aa}', found '{actual_aa}'"
            messages.append(msg)
            raise ValueError(msg)
        
        s[corrected_pos - 1] = new_aa
        messages.append(f"✓ {orig_aa}{orig_pos}{new_aa} (applied at position {corrected_pos})")
    
    return ''.join(s), messages


def retry_mutation_with_pattern(wildtype_seq: str, mutations: List[Tuple[str, int, str]], 
                                cache: UniProtCache = None, uniprotid: str = None, 
                                max_isoforms: int = 6) -> Dict[str, Any]:
    """
    Attempt to retry mutation using different strategies.
    
    Strategy:
    1. Single mutation: Try simple offset search (±10 positions)
    2. Multiple mutations: Use pattern matching
    3. If failed: Try different isoforms (1-6) from UniProt
    
    Args:
        wildtype_seq: Original wildtype sequence
        mutations: List of (fromAA, pos, toAA) tuples
        cache: UniProtCache instance (optional, for isoform fetching)
        uniprotid: UniProt accession for isoform retry
        max_isoforms: Maximum number of isoforms to try
    
    Returns a dict with:
        - success: bool
        - mutant_sequence: str or None
        - offset: int or None
        - messages: List[str]
        - isoform_used: int or None
    """
    result = {
        'success': False,
        'mutant_sequence': None,
        'offset': None,
        'messages': [],
        'isoform_used': None
    }
    
    # Strategy 1: Single mutation - use simple offset search
    if len(mutations) == 1:
        orig_aa, pos, new_aa = mutations[0]
        result['messages'].append(f"Single mutation detected: {orig_aa}{pos}{new_aa}")
        result['messages'].append(f"Using simple offset search (±10 positions)")
        
        success, actual_pos, msg = try_single_mutation_with_offset(
            wildtype_seq, orig_aa, pos, new_aa, max_offset=10
        )
        
        if success:
            # Apply mutation
            s = list(wildtype_seq)
            s[actual_pos - 1] = new_aa
            result['mutant_sequence'] = ''.join(s)
            result['offset'] = actual_pos - pos
            result['success'] = True
            result['messages'].append(f"✓ {msg}")
            return result
        else:
            result['messages'].append(f"✗ {msg}")
            # Continue to isoform retry
    
    # Strategy 2: Multiple mutations - use pattern matching
    else:
        result['messages'].append(f"Multiple mutations detected ({len(mutations)} mutations)")
        result['messages'].append(f"Using pattern matching strategy")
        
        try:
            # Build regex pattern
            pattern, min_pos, max_pos = build_regex_pattern(mutations)
            result['messages'].append(f"Built pattern: {pattern} (positions {min_pos}-{max_pos})")
            
            # Find pattern in sequence
            found_pos = find_pattern_in_sequence(wildtype_seq, pattern, min_pos)
            
            if found_pos:
                # Calculate offset
                offset = found_pos - min_pos
                result['offset'] = offset
                
                if offset == 0:
                    result['messages'].append(f"✓ Pattern found at expected position {found_pos} (no offset)")
                else:
                    result['messages'].append(f"✓ Pattern found at position {found_pos} (offset: {offset:+d})")
                
                # Apply mutations with offset
                mutant_seq, apply_messages = apply_mutations_with_offset(wildtype_seq, mutations, offset)
                result['messages'].extend(apply_messages)
                result['mutant_sequence'] = mutant_seq
                result['success'] = True
                return result
            else:
                result['messages'].append(f"✗ Pattern not found in sequence")
                # Continue to isoform retry
        
        except Exception as e:
            result['messages'].append(f"✗ Pattern matching error: {str(e)}")
            # Continue to isoform retry
    
    # Strategy 3: Try different isoforms (using cache if available)
    if uniprotid and max_isoforms > 0 and cache:
        result['messages'].append(f"\nTrying alternative isoforms...")
        
        # Extract primary UniProt ID
        primary_id = get_primary_uniprot_id(uniprotid)
        
        for iso_num in range(1, max_isoforms + 1):
            try:
                result['messages'].append(f"  Attempting isoform {iso_num}...")
                
                # Fetch isoform sequence using cache
                iso_seq = cache.fetch_isoform_sequence(primary_id, iso_num)
                if not iso_seq:
                    result['messages'].append(f"  ✗ Isoform {iso_num} not found")
                    continue
                
                result['messages'].append(f"  ✓ Fetched isoform {iso_num} (length: {len(iso_seq)})")
                
                # Retry with isoform sequence
                if len(mutations) == 1:
                    # Single mutation with offset search
                    orig_aa, pos, new_aa = mutations[0]
                    success, actual_pos, msg = try_single_mutation_with_offset(
                        iso_seq, orig_aa, pos, new_aa, max_offset=10
                    )
                    
                    if success:
                        s = list(iso_seq)
                        s[actual_pos - 1] = new_aa
                        result['mutant_sequence'] = ''.join(s)
                        result['offset'] = actual_pos - pos
                        result['success'] = True
                        result['isoform_used'] = iso_num
                        result['messages'].append(f"  ✓ SUCCESS with isoform {iso_num}: {msg}")
                        return result
                
                else:
                    # Multiple mutations with pattern matching
                    pattern, min_pos, max_pos = build_regex_pattern(mutations)
                    found_pos = find_pattern_in_sequence(iso_seq, pattern, min_pos)
                    
                    if found_pos:
                        offset = found_pos - min_pos
                        mutant_seq, apply_messages = apply_mutations_with_offset(iso_seq, mutations, offset)
                        
                        result['mutant_sequence'] = mutant_seq
                        result['offset'] = offset
                        result['success'] = True
                        result['isoform_used'] = iso_num
                        result['messages'].append(f"  ✓ SUCCESS with isoform {iso_num} (offset: {offset:+d})")
                        return result
                
                result['messages'].append(f"  ✗ Isoform {iso_num} did not match")
            
            except Exception as e:
                result['messages'].append(f"  ✗ Isoform {iso_num} error: {str(e)}")
                continue
        
        result['messages'].append(f"\n✗ All isoforms (1-{max_isoforms}) failed")
    
    return result


# -------- File processing --------
def list_yaml_paths(root: str) -> List[str]:
    if os.path.isdir(root):
        return sorted(glob.glob(os.path.join(root, "**", "*.yaml"), recursive=True))
    return [root] if root.lower().endswith(".yaml") else []


def cluster_files_by_uniprot(yaml_paths: List[str]) -> Dict[str, List[Tuple[str, Any]]]:
    """Cluster YAML files by UniProt ID."""
    clusters = defaultdict(list)
    
    for path in yaml_paths:
        try:
            with open(path, "r", encoding="utf-8") as fh:
                doc = yaml.safe_load(fh)
            
            if not doc or not isinstance(doc, dict):
                continue
            
            proteins = doc.get("proteins")
            if not isinstance(proteins, list):
                continue
            
            for prot in proteins:
                if not isinstance(prot, dict):
                    continue
                
                # Get the original uniprotid (may contain multiple IDs)
                original_id = get_original_uniprot_id(prot)
                if not original_id:
                    continue
                
                # Extract primary (first) ID for clustering
                primary_id = get_primary_uniprot_id(original_id)
                if not primary_id:
                    continue
                
                # Store the protein with its primary UniProt ID
                clusters[primary_id].append((path, prot))
        
        except Exception as e:
            print(f"[WARNING] Error reading {path}: {e}", file=sys.stderr)
    
    return dict(clusters)


def process_issues_directory(input_dir: str, output_dir: str, report_path: str, 
                            checkpoint: CheckpointManager, cache: UniProtCache,
                            resume: bool = False):
    """
    Process all YAML files in the issues directory and retry mutations.
    
    **NEW**: Automatically fetch missing wildtype_sequence from UniProt
    
    Features:
    - Checkpoint support for resumable processing
    - Organized output by UniProt ID (each in separate folder)
    - Detailed progress output
    - Automatic sequence retrieval with caching
    """
    start_time = time.time()
    
    # Initialize or resume
    if not resume:
        checkpoint.clear()
        print("[INFO] Starting fresh (checkpoint cleared)")
    else:
        stats = checkpoint.get_stats()
        print(f"[INFO] Resuming from checkpoint")
        print(f"       Previously processed: {len(checkpoint.data['processed_files'])} files")
        print(f"       Stats: {stats['success']} success, {stats['failed']} failed, {stats['skipped']} skipped")
    
    # List all YAML files
    print_header("SCANNING FILES")
    yaml_files = list_yaml_paths(input_dir)
    print(f"Found {len(yaml_files)} YAML files")
    
    if not yaml_files:
        print("[WARNING] No YAML files found")
        return
    
    # Cluster files by UniProt ID
    print_header("CLUSTERING BY UNIPROT ID")
    clusters = cluster_files_by_uniprot(yaml_files)
    print(f"Organized into {len(clusters)} UniProt ID clusters")
    
    # Count total entries to process
    total_entries = sum(len(entries) for entries in clusters.values())
    if resume:
        # Subtract already processed entries
        processed_entries = len(checkpoint.data['processed_files'])
        remaining_entries = total_entries - processed_entries
        print(f"Total entries: {total_entries} ({remaining_entries} remaining)")
    else:
        print(f"Total entries: {total_entries}")
    
    # Process clusters
    print_header("PROCESSING MUTATIONS")
    report_rows = []
    processed_count = 0
    
    for cluster_idx, (uniprotid, entries) in enumerate(sorted(clusters.items()), 1):
        # Create output directory for this UniProt ID
        uniprot_output_dir = os.path.join(output_dir, uniprotid)
        os.makedirs(uniprot_output_dir, exist_ok=True)
        
        cluster_success = 0
        cluster_failed = 0
        cluster_skipped = 0
        
        print(f"\n[{cluster_idx}/{len(clusters)}] Processing {uniprotid} ({len(entries)} entries)")
        
        for entry_idx, (file_path, prot) in enumerate(entries, 1):
            # Skip if already processed (for resume)
            if resume and checkpoint.is_processed(file_path):
                continue
            
            # Display progress
            print_progress(entry_idx, len(entries), 
                         prefix=f"    {uniprotid}", 
                         suffix=f"({entry_idx}/{len(entries)})")
            
            # Extract data
            basename = os.path.basename(file_path)
            prot_id = prot.get("id", "")
            variant = prot.get("variant_description", "")
            wildtype_seq = prot.get("wildtype_sequence", "")
            original_uniprotid = get_original_uniprot_id(prot)
            
            # **NEW**: Auto-fetch missing wildtype_sequence
            if not wildtype_seq:
                primary_uniprotid = get_primary_uniprot_id(original_uniprotid)
                
                if primary_uniprotid:
                    print(f"\n    [AUTO-FETCH] Missing wildtype_sequence for {prot_id}")
                    print(f"                 Attempting to fetch from UniProt ({primary_uniprotid})...")
                    
                    # Try canonical first
                    fetched_seq = cache.fetch_canonical_sequence(primary_uniprotid)
                    
                    if fetched_seq:
                        wildtype_seq = fetched_seq
                        checkpoint.update_stats(auto_fetched=1)
                        print(f"                 ✓ Successfully fetched canonical sequence (length: {len(wildtype_seq)})")
                    else:
                        # Try isoforms
                        print(f"                 ✗ Canonical not found, trying isoforms...")
                        for iso_num in range(1, 7):
                            iso_seq = cache.fetch_isoform_sequence(primary_uniprotid, iso_num)
                            if iso_seq:
                                wildtype_seq = iso_seq
                                checkpoint.update_stats(auto_fetched=1)
                                print(f"                 ✓ Successfully fetched isoform {iso_num} (length: {len(wildtype_seq)})")
                                break
            
            # Check if we still don't have a sequence
            if not wildtype_seq:
                checkpoint.mark_processed(file_path, success=False)
                checkpoint.update_stats(skipped=1)
                cluster_skipped += 1
                report_rows.append([
                    file_path, original_uniprotid, prot_id, variant or '', "SKIPPED", 
                    '', '', '', '', "No wildtype_sequence (auto-fetch failed)"
                ])
                continue
            
            if not variant:
                checkpoint.mark_processed(file_path, success=False)
                checkpoint.update_stats(skipped=1)
                cluster_skipped += 1
                report_rows.append([
                    file_path, original_uniprotid, prot_id, '', "SKIPPED", 
                    '', len(wildtype_seq), '', '', "No variant_description"
                ])
                continue
            
            # Parse mutations
            try:
                mutations = parse_variant_description(variant)
                if not mutations:
                    checkpoint.mark_processed(file_path, success=False)
                    checkpoint.update_stats(skipped=1)
                    cluster_skipped += 1
                    report_rows.append([
                        file_path, original_uniprotid, prot_id, variant, "SKIPPED", 
                        '', len(wildtype_seq), '', '', "No mutations parsed"
                    ])
                    continue
            except Exception as e:
                checkpoint.mark_processed(file_path, success=False)
                checkpoint.update_stats(failed=1)
                cluster_failed += 1
                report_rows.append([
                    file_path, original_uniprotid, prot_id, variant, "ERROR", 
                    '', len(wildtype_seq), '', '', f"Parse error: {str(e)}"
                ])
                continue
            
            # Retry mutation with pattern matching and isoform fallback
            try:
                # Get primary UniProt ID for isoform fetching
                primary_uniprotid = get_primary_uniprot_id(original_uniprotid)
                
                result = retry_mutation_with_pattern(
                    wildtype_seq, 
                    mutations, 
                    cache=cache,
                    uniprotid=primary_uniprotid,
                    max_isoforms=6
                )
                
                if result['success']:
                    status = "SUCCESS"
                    checkpoint.update_stats(success=1)
                    checkpoint.mark_processed(file_path, success=True)
                    cluster_success += 1
                    
                    # Update YAML and save to UniProt-specific directory
                    with open(file_path, "r", encoding="utf-8") as fh:
                        doc = yaml.safe_load(fh)
                    
                    # Find and update the protein
                    for i, p in enumerate(doc.get("proteins", [])):
                        if p.get("id") == prot_id:
                            p["mutant_sequence"] = result['mutant_sequence']
                            p["position_offset"] = result['offset']
                            p["retry_status"] = "pattern_matched"
                            # Add isoform info if used
                            if result.get('isoform_used'):
                                p["isoform_used"] = result['isoform_used']
                                p["retry_status"] = f"isoform_{result['isoform_used']}_matched"
                            break
                    
                    # Save updated YAML
                    out_path = os.path.join(uniprot_output_dir, basename)
                    with open(out_path, "w", encoding="utf-8") as fw:
                        yaml.safe_dump(doc, fw, allow_unicode=True, sort_keys=False)
                    
                else:
                    status = "FAILED"
                    checkpoint.update_stats(failed=1)
                    checkpoint.mark_processed(file_path, success=False)
                    cluster_failed += 1
                
                # Include isoform info in report
                isoform_info = f"isoform_{result.get('isoform_used')}" if result.get('isoform_used') else ""
                
                report_rows.append([
                    file_path, original_uniprotid, prot_id, variant, status,
                    result.get('offset', ''), len(wildtype_seq),
                    len(result['mutant_sequence']) if result['mutant_sequence'] else '',
                    isoform_info,
                    '\n'.join(result['messages'])
                ])
                
            except Exception as e:
                status = "ERROR"
                checkpoint.update_stats(failed=1)
                checkpoint.mark_processed(file_path, success=False)
                cluster_failed += 1
                report_rows.append([
                    file_path, original_uniprotid, prot_id, variant or '', "ERROR",
                    '', len(wildtype_seq) if wildtype_seq else '', '', '', str(e)
                ])
            
            processed_count += 1
            
            # Save checkpoint periodically (every 10 files)
            if processed_count % 10 == 0:
                checkpoint.save()
        
        # Print cluster summary
        print(f"\n    ✓ {uniprotid} complete: {cluster_success} success, {cluster_failed} failed, {cluster_skipped} skipped")
    
    # Final checkpoint save
    checkpoint.save()
    
    # Save report
    os.makedirs(os.path.dirname(report_path), exist_ok=True)
    with open(report_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow([
            "file_path", "uniprotid", "protein_id", "variant_description",
            "status", "position_offset", "wildtype_length", "mutant_length", "isoform", "messages"
        ])
        writer.writerows(report_rows)
    
    # Print cache statistics
    cache.print_stats()
    
    # Print final summary
    elapsed_time = time.time() - start_time
    stats = checkpoint.get_stats()
    
    print_header("FINAL SUMMARY")
    print(f"Total processing time: {format_time(elapsed_time)}")
    print(f"\nResults:")
    print(f"  ✓ Success:      {stats['success']:4d} entries")
    print(f"  ✗ Failed:       {stats['failed']:4d} entries")
    print(f"  - Skipped:      {stats['skipped']:4d} entries")
    print(f"  🔄 Auto-fetched: {stats['auto_fetched']:4d} sequences")
    print(f"  ━━━━━━━━━━━━━━━━━━━━━━━━━━━")
    print(f"  = Total:        {total_entries:4d} entries")
    
    if stats['success'] > 0:
        success_rate = 100 * stats['success'] / (stats['success'] + stats['failed'])
        print(f"\nSuccess rate: {success_rate:.1f}%")
    
    if stats['auto_fetched'] > 0:
        print(f"\n🎉 Rescued {stats['auto_fetched']} entries by auto-fetching sequences!")
    
    print(f"\nOutput locations:")
    print(f"  📁 Updated YAMLs: {os.path.abspath(output_dir)}/")
    print(f"     (Organized by UniProt ID)")
    print(f"  📄 Report:        {os.path.abspath(report_path)}")
    print(f"  💾 Checkpoint:    {os.path.abspath(checkpoint.checkpoint_file)}")
    print()


# -------- Main --------
def main():
    ap = argparse.ArgumentParser(
        description="Enhanced mutation retry with automatic sequence fetching",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Start fresh
  python step1_enhanced.py --input ./mutation_exports/issues
  
  # Resume from checkpoint (after interruption)
  python step1_enhanced.py --input ./mutation_exports/issues --resume
  
  # Custom cache directory
  python step1_enhanced.py --input ./mutation_exports/issues --cache-dir ./my_cache
  
  # Clear checkpoint and start over
  python step1_enhanced.py --input ./mutation_exports/issues --clear-checkpoint

Features:
  - Auto-fetches missing wildtype_sequence from UniProt
  - Tries canonical sequence first, then isoforms (1-6)
  - Local disk cache to avoid redundant API calls
  - Rate limiting with exponential backoff retry
  - Checkpoint support for resumable processing
        """
    )
    ap.add_argument("--input", "-i", required=True, 
                    help="Issues directory containing YAML files")
    ap.add_argument("--output", "-o", default="./mutation_exports/retry_success",
                    help="Output directory for corrected YAMLs (organized by UniProt ID)")
    ap.add_argument("--report", "-r", default="./mutation_reports/retry_report.csv",
                    help="Path for retry report CSV")
    ap.add_argument("--checkpoint", "-c", default="./mutation_reports/.checkpoint.json",
                    help="Path for checkpoint file (for resumable processing)")
    ap.add_argument("--cache-dir", default="./uniprot_cache",
                    help="Directory for UniProt sequence cache")
    ap.add_argument("--rate-limit", type=float, default=1.0,
                    help="Minimum seconds between UniProt API requests (default: 1.0)")
    ap.add_argument("--resume", action="store_true",
                    help="Resume from last checkpoint")
    ap.add_argument("--clear-checkpoint", action="store_true",
                    help="Clear checkpoint and start fresh")
    args = ap.parse_args()
    
    # Convert to absolute paths
    args.input = os.path.abspath(args.input)
    args.output = os.path.abspath(args.output)
    args.report = os.path.abspath(args.report)
    args.checkpoint = os.path.abspath(args.checkpoint)
    args.cache_dir = os.path.abspath(args.cache_dir)
    
    if not os.path.exists(args.input):
        print(f"[ERROR] Input directory does not exist: {args.input}", file=sys.stderr)
        sys.exit(1)
    
    # Initialize checkpoint manager
    checkpoint = CheckpointManager(args.checkpoint)
    
    # Clear checkpoint if requested
    if args.clear_checkpoint:
        checkpoint.clear()
        print("[INFO] Checkpoint cleared")
        return
    
    # Initialize UniProt cache
    cache = UniProtCache(cache_dir=args.cache_dir, rate_limit=args.rate_limit)
    
    # Print banner
    print("=" * 70)
    print("  ENHANCED MUTATION RETRY WITH AUTO SEQUENCE FETCHING")
    print("=" * 70)
    print(f"\nConfiguration:")
    print(f"  Input:      {args.input}")
    print(f"  Output:     {args.output}")
    print(f"  Report:     {args.report}")
    print(f"  Checkpoint: {args.checkpoint}")
    print(f"  Cache dir:  {args.cache_dir}")
    print(f"  Rate limit: {args.rate_limit}s between API requests")
    print(f"  Mode:       {'RESUME' if args.resume else 'FRESH START'}")
    
    # Process
    process_issues_directory(args.input, args.output, args.report, checkpoint, cache, args.resume)


if __name__ == "__main__":
    main()
