#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Retry failed mutations using pattern-based sequence alignment.

This script processes YAML files from the issues directory and attempts to:
1. Cluster files by UniProt ID (organized into separate folders)
2. Build regex patterns from mutation positions to locate the correct region in the sequence
3. Calculate position offset and re-apply mutations with corrected positions

Features:
- Checkpoint support: resume from last interruption
- Detailed progress output
- Organized output by UniProt ID

Usage:
    python retry_mutations_with_pattern.py --input ./mutation_exports/issues
    
    # Resume from checkpoint
    python retry_mutations_with_pattern.py --input ./mutation_exports/issues --resume
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
from typing import Dict, List, Tuple, Optional, Any

import yaml


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
                'skipped': 0
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
            'stats': {'total': 0, 'success': 0, 'failed': 0, 'skipped': 0},
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
    """
    Check if a uniprotid string contains multiple IDs.
    
    Args:
        uniprotid: String that may contain one or more UniProt IDs
    
    Returns:
        True if multiple IDs are detected, False otherwise
    """
    if not uniprotid:
        return False
    
    uniprotid = str(uniprotid).strip()
    
    # Check for common separators
    for sep in [';', ',', '|']:
        if sep in uniprotid:
            return True
    
    return False


def get_primary_uniprot_id(uniprotid: str) -> str:
    """
    Extract the primary (first) UniProt ID from a string that may contain multiple IDs.
    
    Args:
        uniprotid: String that may contain one or more UniProt IDs separated by ;, , or |
    
    Returns:
        The first/primary UniProt ID
    """
    if not uniprotid:
        return ""
    
    uniprotid = str(uniprotid).strip()
    
    # Split by common separators and take the first ID
    for sep in [';', ',', '|']:
        if sep in uniprotid:
            return uniprotid.split(sep)[0].strip()
    
    return uniprotid


def get_original_uniprot_id(prot: Dict) -> str:
    """
    Get the original uniprotid from protein dict (may contain multiple IDs).
    
    Returns:
        The original uniprotid string as stored in the YAML
    """
    return str(prot.get("uniprotid") or prot.get("uniprot_id") or prot.get("uniprot") or "")


# -------- Parse variant_description (same as original) --------
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
def fetch_uniprot_isoform(accession: str, isoform: int, timeout: int = 30) -> str:
    """
    Fetch a specific isoform sequence from UniProt.
    
    Args:
        accession: UniProt accession (e.g., 'P77256')
        isoform: Isoform number (1, 2, 3, etc.)
        timeout: Request timeout in seconds
    
    Returns:
        Amino acid sequence string
    
    Raises:
        Exception if fetch fails
    """
    import requests
    
    query = f"accession:{accession}-{isoform}"
    url = f"https://rest.uniprot.org/uniprotkb/stream?format=fasta&query={query}"
    
    try:
        r = requests.get(url, timeout=timeout)
        r.raise_for_status()
        text = r.text.strip()
        
        if not text.startswith(">"):
            raise ValueError(f"Invalid FASTA for {accession}-{isoform}")
        
        seq = "".join(line.strip() for line in text.splitlines()[1:] if line and not line.startswith(">"))
        
        if not seq:
            raise ValueError(f"No sequence parsed for {accession}-{isoform}")
        
        return seq
    except Exception as e:
        raise Exception(f"Failed to fetch isoform {isoform}: {str(e)}")


def try_single_mutation_with_offset(seq: str, orig_aa: str, pos: int, new_aa: str, 
                                     max_offset: int = 10) -> Tuple[bool, int, str]:
    """
    Try to apply a single mutation with position offset search.
    
    Args:
        seq: Sequence to mutate
        orig_aa: Original amino acid expected
        pos: Original position (1-based)
        new_aa: New amino acid
        max_offset: Maximum offset to search (± offset)
    
    Returns:
        (success, actual_position, message)
    """
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
    """
    Build a regex pattern from mutations to locate the region in the sequence.
    
    Returns:
        - pattern: regex pattern string
        - min_pos: minimum position in mutations
        - max_pos: maximum position in mutations
    """
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
    """
    Find the pattern in the sequence, searching around the expected position.
    
    Returns:
        - The 1-based starting position if found, None otherwise
    """
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
    """
    Apply mutations with position offset.
    
    Returns:
        - mutated sequence
        - list of verification messages
    """
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
                                uniprotid: str = None, max_isoforms: int = 6) -> Dict[str, Any]:
    """
    Attempt to retry mutation using different strategies.
    
    Strategy:
    1. Single mutation: Try simple offset search (±10 positions)
    2. Multiple mutations: Use pattern matching
    3. If failed: Try different isoforms (1-6) from UniProt
    
    Args:
        wildtype_seq: Original wildtype sequence
        mutations: List of (fromAA, pos, toAA) tuples
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
    
    # Strategy 3: Try different isoforms
    if uniprotid and max_isoforms > 0:
        result['messages'].append(f"\nTrying alternative isoforms...")
        
        # Extract primary UniProt ID
        primary_id = get_primary_uniprot_id(uniprotid)
        
        import time
        
        for iso_num in range(1, max_isoforms + 1):
            try:
                result['messages'].append(f"  Attempting isoform {iso_num}...")
                
                # Fetch isoform sequence
                iso_seq = fetch_uniprot_isoform(primary_id, iso_num, timeout=30)
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
                
                # Add delay between isoform requests to be polite to UniProt API
                time.sleep(0.2)
            
            except Exception as e:
                result['messages'].append(f"  ✗ Isoform {iso_num} error: {str(e)}")
                time.sleep(0.1)  # Brief delay even on error
                continue
        
        result['messages'].append(f"\n✗ All isoforms (1-{max_isoforms}) failed")
    
    return result


# -------- File processing --------
def list_yaml_paths(root: str) -> List[str]:
    if os.path.isdir(root):
        return sorted(glob.glob(os.path.join(root, "**", "*.yaml"), recursive=True))
    return [root] if root.lower().endswith(".yaml") else []


def cluster_files_by_uniprot(yaml_paths: List[str]) -> Dict[str, List[Tuple[str, Any]]]:
    """
    Cluster YAML files by UniProt ID.
    
    If a protein has multiple UniProt IDs (separated by semicolon/comma),
    only the first (primary) ID is used for clustering.
    
    Returns:
        Dict[primary_uniprotid, List[(file_path, protein_data)]]
    """
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
                            checkpoint: CheckpointManager, resume: bool = False):
    """
    Process all YAML files in the issues directory and retry mutations.
    
    Features:
    - Checkpoint support for resumable processing
    - Organized output by UniProt ID (each in separate folder)
    - Detailed progress output
    """
    start_time = time.time()
    
    # Initialize or resume
    if not resume:
        checkpoint.clear()
        print("[INFO] Starting fresh (checkpoint cleared)")
    else:
        stats = checkpoint.get_stats()
        if stats['total'] > 0:
            print(f"[INFO] Resuming from checkpoint:")
            print(f"  - Previously processed: {len(checkpoint.data['processed_files'])} files")
            print(f"  - Success: {stats['success']}, Failed: {stats['failed']}, Skipped: {stats['skipped']}")
    
    # Get all YAML files
    yaml_paths = list_yaml_paths(input_dir)
    if not yaml_paths:
        print("[ERROR] No YAML files found in issues directory.")
        return
    
    print(f"\n[INFO] Found {len(yaml_paths)} YAML files in total")
    
    # Step 1: Cluster by UniProt ID
    print_header("STEP 1: Clustering by UniProt ID")
    
    clusters = cluster_files_by_uniprot(yaml_paths)
    print(f"✓ Clustered into {len(clusters)} unique UniProt IDs")
    print(f"  (Note: When multiple IDs exist, only the first/primary ID is used)\n")
    
    # Show top 10 UniProt IDs by file count
    sorted_clusters = sorted([(uid, len(files)) for uid, files in clusters.items()], 
                            key=lambda x: -x[1])
    print("Top 10 UniProt IDs by file count:")
    for i, (uniprotid, count) in enumerate(sorted_clusters[:10], 1):
        print(f"  {i:2d}. {uniprotid}: {count} files")
    
    # Calculate total work
    total_entries = sum(len(files) for files in clusters.values())
    print(f"\nTotal entries to process: {total_entries}")
    
    # Step 2 & 3: Process each cluster
    print_header("STEP 2 & 3: Pattern Matching and Mutation Retry")
    
    report_rows = []
    os.makedirs(output_dir, exist_ok=True)
    
    # Process counters
    processed_count = 0
    cluster_num = 0
    
    for uniprotid, file_prots in clusters.items():
        cluster_num += 1
        cluster_success = 0
        cluster_failed = 0
        cluster_skipped = 0
        
        # Create UniProt-specific output directory
        uniprot_output_dir = os.path.join(output_dir, uniprotid)
        os.makedirs(uniprot_output_dir, exist_ok=True)
        
        print(f"\n[{cluster_num}/{len(clusters)}] Processing UniProt: {uniprotid} ({len(file_prots)} entries)")
        print(f"    Output directory: {os.path.relpath(uniprot_output_dir)}")
        
        for entry_idx, (file_path, prot) in enumerate(file_prots, 1):
            # Skip if already processed (checkpoint)
            if checkpoint.is_processed(file_path):
                processed_count += 1
                cluster_skipped += 1
                continue
            
            basename = os.path.basename(file_path)
            prot_id = prot.get("id", "unknown")
            variant = prot.get("variant_description") or prot.get("variant") or prot.get("mutation")
            wildtype_seq = prot.get("wildtype_sequence")
            
            # Get original uniprotid (may contain multiple IDs) for reporting
            original_uniprotid = get_original_uniprot_id(prot)
            
            # Skip entries with multiple UniProt IDs
            if has_multiple_uniprot_ids(original_uniprotid):
                checkpoint.mark_processed(file_path, success=False)
                checkpoint.update_stats(skipped=1)
                cluster_skipped += 1
                processed_count += 1
                
                # Show progress
                print_progress(entry_idx, len(file_prots), 
                              prefix=f"    [{uniprotid}]",
                              suffix=f"({entry_idx}/{len(file_prots)}) {basename[:40]}...")
                
                report_rows.append([
                    file_path, original_uniprotid, prot_id, variant or '', "SKIPPED", 
                    '', '', '', '', "Multiple UniProt IDs detected - ambiguous sequence source"
                ])
                continue
            
            # Progress bar
            print_progress(entry_idx, len(file_prots), 
                          prefix=f"    [{uniprotid}]",
                          suffix=f"({entry_idx}/{len(file_prots)}) {basename[:40]}...")
            
            # Skip conditions
            if not wildtype_seq:
                checkpoint.mark_processed(file_path, success=False)
                checkpoint.update_stats(skipped=1)
                cluster_skipped += 1
                report_rows.append([
                    file_path, original_uniprotid, prot_id, variant or '', "SKIPPED", 
                    '', '', '', '', "No wildtype_sequence"
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
    
    # Print final summary
    elapsed_time = time.time() - start_time
    stats = checkpoint.get_stats()
    
    print_header("FINAL SUMMARY")
    print(f"Total processing time: {format_time(elapsed_time)}")
    print(f"\nResults:")
    print(f"  ✓ Success:  {stats['success']:4d} entries")
    print(f"  ✗ Failed:   {stats['failed']:4d} entries")
    print(f"  - Skipped:  {stats['skipped']:4d} entries")
    print(f"  ━━━━━━━━━━━━━━━━━━━━━━")
    print(f"  = Total:    {total_entries:4d} entries")
    
    if stats['success'] > 0:
        success_rate = 100 * stats['success'] / (stats['success'] + stats['failed'])
        print(f"\nSuccess rate: {success_rate:.1f}%")
    
    print(f"\nOutput locations:")
    print(f"  📁 Updated YAMLs: {os.path.abspath(output_dir)}/")
    print(f"     (Organized by UniProt ID)")
    print(f"  📄 Report:        {os.path.abspath(report_path)}")
    print(f"  💾 Checkpoint:    {os.path.abspath(checkpoint.checkpoint_file)}")
    print()


# -------- Main --------
def main():
    ap = argparse.ArgumentParser(
        description="Retry failed mutations using pattern matching with checkpoint support",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Start fresh
  python retry_mutations_with_pattern.py --input ./mutation_exports/issues
  
  # Resume from checkpoint (after interruption)
  python retry_mutations_with_pattern.py --input ./mutation_exports/issues --resume
  
  # Clear checkpoint and start over
  python retry_mutations_with_pattern.py --input ./mutation_exports/issues --clear-checkpoint
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
    
    # Print banner
    print("=" * 70)
    print("  MUTATION RETRY WITH PATTERN MATCHING")
    print("=" * 70)
    print(f"\nConfiguration:")
    print(f"  Input:      {args.input}")
    print(f"  Output:     {args.output}")
    print(f"  Report:     {args.report}")
    print(f"  Checkpoint: {args.checkpoint}")
    print(f"  Mode:       {'RESUME' if args.resume else 'FRESH START'}")
    
    # Process
    process_issues_directory(args.input, args.output, args.report, checkpoint, args.resume)


if __name__ == "__main__":
    main()
