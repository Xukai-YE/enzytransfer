#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Step 2: Multi-UniProt-ID Resolver

This script handles YAML files with multiple UniProt IDs (separated by ;, ,, or |).

Strategy:
1. Identify files with multiple UniProt IDs
2. Fetch sequences for ALL UniProt IDs (canonical + isoforms)
3. Try mutation on each sequence in order
4. If direct match found:
   - Update YAML with the matching UniProt ID
   - Save to resolved directory
5. If no match found:
   - Save to unresolved directory for manual review

Usage:
    python step2_multi.py --input ./mutation_exports/issues
    
    # Resume from checkpoint
    python step2_multi_.py --input ./mutation_exports/issues --resume
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

# Import cache from step1_enhanced
sys.path.insert(0, os.path.dirname(__file__))
from step1 import (
    UniProtCache,
    parse_variant_description,
    get_primary_uniprot_id,
    get_original_uniprot_id,
    has_multiple_uniprot_ids,
    CheckpointManager,
    print_header,
    print_progress,
    format_time
)


# ==================== Multi-ID Specific Functions ====================

def split_uniprot_ids(uniprotid_str: str) -> List[str]:
    """
    Split a string containing multiple UniProt IDs.
    
    Args:
        uniprotid_str: String like "P12345;Q67890" or "A1B2C3,D4E5F6"
    
    Returns:
        List of individual UniProt IDs
    """
    if not uniprotid_str:
        return []
    
    # Try different separators
    for sep in [';', ',', '|']:
        if sep in uniprotid_str:
            ids = [uid.strip() for uid in uniprotid_str.split(sep)]
            return [uid for uid in ids if uid]
    
    # Single ID
    return [uniprotid_str.strip()]


def try_direct_mutation_match(seq: str, mutations: List[Tuple[str, int, str]]) -> Tuple[bool, Optional[str], List[str]]:
    """
    Try to directly apply mutations to a sequence (no offset search).
    
    This is a strict check - mutations must match exactly at specified positions.
    
    Args:
        seq: Sequence to try
        mutations: List of (fromAA, pos, toAA) tuples
    
    Returns:
        (success, mutant_sequence, messages)
    """
    messages = []
    
    # Check all mutations can be applied
    for orig_aa, pos, new_aa in mutations:
        # Check position is in range
        if not (1 <= pos <= len(seq)):
            messages.append(f"✗ Position {pos} out of range (sequence length: {len(seq)})")
            return (False, None, messages)
        
        # Check amino acid matches
        actual_aa = seq[pos - 1]
        if actual_aa != orig_aa:
            messages.append(f"✗ Mismatch at position {pos}: expected '{orig_aa}', found '{actual_aa}'")
            return (False, None, messages)
    
    # All checks passed - apply mutations
    s = list(seq)
    for orig_aa, pos, new_aa in mutations:
        s[pos - 1] = new_aa
        messages.append(f"✓ Applied {orig_aa}{pos}{new_aa}")
    
    mutant_seq = ''.join(s)
    messages.append(f"✓ SUCCESS: All {len(mutations)} mutations applied")
    
    return (True, mutant_seq, messages)


def try_all_uniprot_ids(uniprot_ids: List[str], mutations: List[Tuple[str, int, str]], 
                       cache: UniProtCache, max_isoforms: int = 6) -> Dict[str, Any]:
    """
    Try mutations on all UniProt IDs and their isoforms.
    
    Args:
        uniprot_ids: List of UniProt IDs to try
        mutations: Mutations to apply
        cache: UniProtCache instance
        max_isoforms: Maximum number of isoforms to try per ID
    
    Returns:
        Dict with:
        - success: bool
        - matched_id: str (UniProt ID that worked)
        - matched_isoform: int or None
        - mutant_sequence: str
        - wildtype_sequence: str
        - messages: List[str]
    """
    result = {
        'success': False,
        'matched_id': None,
        'matched_isoform': None,
        'mutant_sequence': None,
        'wildtype_sequence': None,
        'messages': []
    }
    
    result['messages'].append(f"Trying {len(uniprot_ids)} UniProt IDs: {', '.join(uniprot_ids)}")
    
    # Try each UniProt ID in order
    for uid in uniprot_ids:
        result['messages'].append(f"\n--- Testing {uid} ---")
        
        # Try canonical first
        result['messages'].append(f"  Fetching canonical sequence...")
        canonical_seq = cache.fetch_canonical_sequence(uid)
        
        if canonical_seq:
            result['messages'].append(f"  ✓ Got canonical (length: {len(canonical_seq)})")
            result['messages'].append(f"  Trying mutations on canonical...")
            
            success, mutant_seq, mut_messages = try_direct_mutation_match(canonical_seq, mutations)
            result['messages'].extend([f"    {msg}" for msg in mut_messages])
            
            if success:
                result['success'] = True
                result['matched_id'] = uid
                result['matched_isoform'] = None
                result['mutant_sequence'] = mutant_seq
                result['wildtype_sequence'] = canonical_seq
                result['messages'].append(f"\n🎉 MATCH FOUND: {uid} (canonical)")
                return result
        else:
            result['messages'].append(f"  ✗ Canonical not found")
        
        # Try isoforms
        result['messages'].append(f"  Trying isoforms 1-{max_isoforms}...")
        for iso_num in range(1, max_isoforms + 1):
            iso_seq = cache.fetch_isoform_sequence(uid, iso_num)
            
            if not iso_seq:
                continue
            
            result['messages'].append(f"  ✓ Got isoform-{iso_num} (length: {len(iso_seq)})")
            result['messages'].append(f"  Trying mutations on isoform-{iso_num}...")
            
            success, mutant_seq, mut_messages = try_direct_mutation_match(iso_seq, mutations)
            result['messages'].extend([f"    {msg}" for msg in mut_messages])
            
            if success:
                result['success'] = True
                result['matched_id'] = uid
                result['matched_isoform'] = iso_num
                result['mutant_sequence'] = mutant_seq
                result['wildtype_sequence'] = iso_seq
                result['messages'].append(f"\n🎉 MATCH FOUND: {uid} (isoform-{iso_num})")
                return result
    
    result['messages'].append(f"\n✗ No match found in any of the {len(uniprot_ids)} UniProt IDs")
    return result


# ==================== File Processing ====================

def identify_multi_id_files(yaml_paths: List[str]) -> Tuple[List[str], List[str]]:
    """
    Separate YAML files into single-ID and multi-ID categories.
    
    Returns:
        (single_id_files, multi_id_files)
    """
    single_id_files = []
    multi_id_files = []
    
    for path in yaml_paths:
        try:
            with open(path, "r", encoding="utf-8") as fh:
                doc = yaml.safe_load(fh)
            
            if not doc or not isinstance(doc, dict):
                continue
            
            proteins = doc.get("proteins")
            if not isinstance(proteins, list):
                continue
            
            # Check if any protein has multiple IDs
            has_multi_id = False
            for prot in proteins:
                if not isinstance(prot, dict):
                    continue
                
                original_id = get_original_uniprot_id(prot)
                if has_multiple_uniprot_ids(original_id):
                    has_multi_id = True
                    break
            
            if has_multi_id:
                multi_id_files.append(path)
            else:
                single_id_files.append(path)
        
        except Exception as e:
            print(f"[WARNING] Error reading {path}: {e}", file=sys.stderr)
    
    return single_id_files, multi_id_files


def process_multi_id_files(input_dir: str, resolved_dir: str, unresolved_dir: str,
                          report_path: str, checkpoint: CheckpointManager,
                          cache: UniProtCache, resume: bool = False):
    """
    Process YAML files with multiple UniProt IDs.
    
    Strategy:
    1. Identify files with multiple UniProt IDs
    2. For each file:
       - Get all UniProt IDs
       - Fetch sequences for all IDs (canonical + isoforms)
       - Try mutations on each sequence in order
       - If match found: update YAML and save to resolved
       - If no match: save to unresolved
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
        print(f"       Stats: {stats['success']} resolved, {stats['failed']} unresolved")
    
    # List all YAML files
    print_header("SCANNING FILES")
    yaml_files = []
    if os.path.isdir(input_dir):
        yaml_files = sorted(glob.glob(os.path.join(input_dir, "**", "*.yaml"), recursive=True))
    print(f"Found {len(yaml_files)} YAML files")
    
    if not yaml_files:
        print("[WARNING] No YAML files found")
        return
    
    # Identify multi-ID files
    print_header("IDENTIFYING MULTI-ID FILES")
    single_id_files, multi_id_files = identify_multi_id_files(yaml_files)
    
    print(f"Single UniProt ID files:   {len(single_id_files)}")
    print(f"Multiple UniProt ID files: {len(multi_id_files)}")
    
    if not multi_id_files:
        print("\n[INFO] No multi-ID files found. Nothing to process.")
        return
    
    # Count total entries
    total_entries = 0
    for path in multi_id_files:
        try:
            with open(path, "r", encoding="utf-8") as fh:
                doc = yaml.safe_load(fh)
            proteins = doc.get("proteins", [])
            for prot in proteins:
                if has_multiple_uniprot_ids(get_original_uniprot_id(prot)):
                    total_entries += 1
        except Exception:
            pass
    
    print(f"Total multi-ID entries to process: {total_entries}")
    
    # Process files
    print_header("PROCESSING MULTI-ID FILES")
    report_rows = []
    processed_count = 0
    resolved_count = 0
    unresolved_count = 0
    
    os.makedirs(resolved_dir, exist_ok=True)
    os.makedirs(unresolved_dir, exist_ok=True)
    
    for file_idx, file_path in enumerate(multi_id_files, 1):
        # Skip if already processed
        if resume and checkpoint.is_processed(file_path):
            continue
        
        # Display progress
        print_progress(file_idx, len(multi_id_files),
                     prefix="Processing",
                     suffix=f"({file_idx}/{len(multi_id_files)})")
        
        try:
            # Read YAML
            with open(file_path, "r", encoding="utf-8") as fh:
                doc = yaml.safe_load(fh)
            
            if not doc or not isinstance(doc, dict):
                continue
            
            proteins = doc.get("proteins", [])
            if not proteins:
                continue
            
            # Track if any protein in this file was resolved
            file_resolved = False
            file_modified = False
            
            # Process each protein
            for prot in proteins:
                if not isinstance(prot, dict):
                    continue
                
                # Get UniProt IDs
                original_uniprotid = get_original_uniprot_id(prot)
                if not has_multiple_uniprot_ids(original_uniprotid):
                    continue  # Skip single-ID proteins
                
                prot_id = prot.get("id", "")
                variant = prot.get("variant_description", "")
                
                print(f"\n\n  Processing protein: {prot_id}")
                print(f"  Original IDs: {original_uniprotid}")
                
                # Parse mutations
                if not variant:
                    print(f"  ✗ No variant_description, skipping")
                    report_rows.append([
                        file_path, original_uniprotid, prot_id, "", "SKIPPED",
                        "", "", "", "No variant_description"
                    ])
                    continue
                
                try:
                    mutations = parse_variant_description(variant)
                    if not mutations:
                        print(f"  ✗ No mutations parsed from '{variant}'")
                        report_rows.append([
                            file_path, original_uniprotid, prot_id, variant, "SKIPPED",
                            "", "", "", "No mutations parsed"
                        ])
                        continue
                    
                    print(f"  Parsed {len(mutations)} mutation(s): {variant}")
                
                except Exception as e:
                    print(f"  ✗ Error parsing mutations: {e}")
                    report_rows.append([
                        file_path, original_uniprotid, prot_id, variant, "ERROR",
                        "", "", "", f"Parse error: {str(e)}"
                    ])
                    continue
                
                # Split into individual UniProt IDs
                uniprot_ids = split_uniprot_ids(original_uniprotid)
                print(f"  Split into {len(uniprot_ids)} IDs: {uniprot_ids}")
                
                # Try all UniProt IDs
                result = try_all_uniprot_ids(uniprot_ids, mutations, cache, max_isoforms=6)
                
                # Print result messages (indented)
                for msg in result['messages']:
                    print(f"  {msg}")
                
                if result['success']:
                    print(f"\n  ✓ RESOLVED: Using {result['matched_id']}")
                    
                    # Update protein in YAML
                    prot['uniprotid'] = result['matched_id']  # Update to matched ID only
                    prot['wildtype_sequence'] = result['wildtype_sequence']
                    prot['mutant_sequence'] = result['mutant_sequence']
                    prot['resolution_status'] = 'multi_id_resolved'
                    if result['matched_isoform']:
                        prot['matched_isoform'] = result['matched_isoform']
                    
                    file_resolved = True
                    file_modified = True
                    resolved_count += 1
                    
                    # Report row
                    isoform_info = f"isoform-{result['matched_isoform']}" if result['matched_isoform'] else "canonical"
                    report_rows.append([
                        file_path, original_uniprotid, prot_id, variant, "RESOLVED",
                        result['matched_id'], isoform_info,
                        len(result['wildtype_sequence']),
                        '\n'.join(result['messages'][:10])  # First 10 lines
                    ])
                
                else:
                    print(f"\n  ✗ UNRESOLVED: No match found")
                    unresolved_count += 1
                    
                    # Report row
                    report_rows.append([
                        file_path, original_uniprotid, prot_id, variant, "UNRESOLVED",
                        "", "", "",
                        '\n'.join(result['messages'][:10])  # First 10 lines
                    ])
                
                processed_count += 1
            
            # Save YAML to appropriate directory
            basename = os.path.basename(file_path)
            
            if file_resolved and file_modified:
                # Save to resolved directory
                out_path = os.path.join(resolved_dir, basename)
                with open(out_path, "w", encoding="utf-8") as fw:
                    yaml.safe_dump(doc, fw, allow_unicode=True, sort_keys=False)
                checkpoint.update_stats(success=1)
                checkpoint.mark_processed(file_path, success=True)
            else:
                # Save to unresolved directory
                out_path = os.path.join(unresolved_dir, basename)
                # Copy original file
                import shutil
                shutil.copy2(file_path, out_path)
                checkpoint.update_stats(failed=1)
                checkpoint.mark_processed(file_path, success=False)
            
            # Save checkpoint periodically
            if processed_count % 10 == 0:
                checkpoint.save()
        
        except Exception as e:
            print(f"\n[ERROR] Failed to process {file_path}: {e}")
            import traceback
            traceback.print_exc()
            checkpoint.mark_processed(file_path, success=False)
            checkpoint.update_stats(failed=1)
    
    # Final checkpoint save
    checkpoint.save()
    
    # Save report
    print("\n")
    print_header("SAVING REPORT")
    os.makedirs(os.path.dirname(report_path), exist_ok=True)
    with open(report_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow([
            "file_path", "original_uniprotid", "protein_id", "variant_description",
            "status", "matched_id", "matched_isoform", "sequence_length", "messages"
        ])
        writer.writerows(report_rows)
    
    print(f"Report saved: {report_path}")
    
    # Print cache statistics
    cache.print_stats()
    
    # Print final summary
    elapsed_time = time.time() - start_time
    stats = checkpoint.get_stats()
    
    print_header("FINAL SUMMARY")
    print(f"Total processing time: {format_time(elapsed_time)}")
    print(f"\nResults:")
    print(f"  ✓ Resolved:   {resolved_count:4d} entries")
    print(f"  ✗ Unresolved: {unresolved_count:4d} entries")
    print(f"  ━━━━━━━━━━━━━━━━━━━━━━━━━━━")
    print(f"  = Total:      {processed_count:4d} entries")
    
    if processed_count > 0:
        resolution_rate = 100 * resolved_count / processed_count
        print(f"\nResolution rate: {resolution_rate:.1f}%")
    
    print(f"\nOutput locations:")
    print(f"  ✅ Resolved:   {os.path.abspath(resolved_dir)}/")
    print(f"  ❌ Unresolved: {os.path.abspath(unresolved_dir)}/")
    print(f"  📄 Report:     {os.path.abspath(report_path)}")
    print(f"  💾 Checkpoint: {os.path.abspath(checkpoint.checkpoint_file)}")
    print()


# ==================== Main ====================

def main():
    ap = argparse.ArgumentParser(
        description="Resolve multi-UniProt-ID entries by trying mutations on all IDs",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Process multi-ID files
  python step2_multi_id_resolver.py --input ./mutation_exports/issues
  
  # Resume from checkpoint
  python step2_multi_id_resolver.py --input ./mutation_exports/issues --resume
  
  # Custom output directories
  python step2_multi_id_resolver.py \
    --input ./mutation_exports/issues \
    --resolved ./mutation_exports/multi_id_resolved \
    --unresolved ./mutation_exports/multi_id_unresolved

How it works:
  1. Identifies YAML files with multiple UniProt IDs (separated by ; , |)
  2. Fetches sequences for ALL IDs (canonical + isoforms 1-6)
  3. Tries mutations on each sequence in order
  4. If direct match found:
     - Updates YAML with the matching UniProt ID
     - Saves to resolved directory
  5. If no match found:
     - Saves to unresolved directory for manual review
        """
    )
    ap.add_argument("--input", "-i", required=True,
                    help="Input directory containing YAML files")
    ap.add_argument("--resolved", default="./mutation_exports/multi_id_resolved",
                    help="Output directory for resolved files")
    ap.add_argument("--unresolved", default="./mutation_exports/multi_id_unresolved",
                    help="Output directory for unresolved files")
    ap.add_argument("--report", "-r", default="./mutation_reports/multi_id_report.csv",
                    help="Path for report CSV")
    ap.add_argument("--checkpoint", "-c", default="./mutation_reports/.checkpoint_multi_id.json",
                    help="Path for checkpoint file")
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
    args.resolved = os.path.abspath(args.resolved)
    args.unresolved = os.path.abspath(args.unresolved)
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
    print("  MULTI-UNIPROT-ID RESOLVER")
    print("=" * 70)
    print(f"\nConfiguration:")
    print(f"  Input:      {args.input}")
    print(f"  Resolved:   {args.resolved}")
    print(f"  Unresolved: {args.unresolved}")
    print(f"  Report:     {args.report}")
    print(f"  Checkpoint: {args.checkpoint}")
    print(f"  Cache dir:  {args.cache_dir}")
    print(f"  Rate limit: {args.rate_limit}s between API requests")
    print(f"  Mode:       {'RESUME' if args.resume else 'FRESH START'}")
    
    # Process
    process_multi_id_files(args.input, args.resolved, args.unresolved,
                          args.report, checkpoint, cache, args.resume)


if __name__ == "__main__":
    main()
