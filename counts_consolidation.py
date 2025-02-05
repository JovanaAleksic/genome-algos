import pandas as pd
import numpy as np
import os
from typing import List, Dict, Set, Tuple
from collections import defaultdict
from pathlib import Path

def fragments(filename):
    fragments = []
    with open(filename, 'r') as f:
        next(f)  # Skip header line 'fragment_id'
        for line in f:
            fragments.append(line.strip())
    return fragments

def build_fragment_lookup(real_fragments_list: List[str], max_offset: int = 2) -> Dict[str, List[Tuple[str, int]]]:
    """Build an efficient lookup dictionary for fragments."""
    lookup = defaultdict(list)
    for frag in real_fragments_list:
        protein, pos = frag.split(':')
        pos = int(pos)
        # Create range of positions for this protein
        for offset in range(-max_offset, max_offset + 1):
            lookup[protein].append((protein, pos + offset))
    return lookup

def find_matching_fragment(protein: str, pos: int, fragment_lookup: Dict[str, List[Tuple[str, int]]]) -> Tuple[str, int]:
    """Find matching real fragment using the lookup dictionary."""
    if protein not in fragment_lookup:
        return None
    
    candidates = fragment_lookup[protein]
    for candidate_protein, candidate_pos in candidates:
        if candidate_pos == pos:
            # Return the base position (without offset)
            base_pos = None
            for offset in range(3):  # max_offset + 1
                if (candidate_protein, candidate_pos - offset) in candidates:
                    base_pos = candidate_pos - offset
                    break
                if (candidate_protein, candidate_pos + offset) in candidates:
                    base_pos = candidate_pos + offset
                    break
            if base_pos is not None:
                return (candidate_protein, base_pos)
    return None

def process_fragment_pairs(df: pd.DataFrame, fragment_lookup: Dict[str, List[Tuple[str, int]]]) -> Dict[str, str]:
    """Process all fragment pairs to create a mapping to real pairs."""
    mapping = {}
    
    # Extract unique fragment pairs for processing
    unique_pairs = df['ID'].unique()
    
    for pair in unique_pairs:
        parts = pair.split(':')
        if len(parts) != 4:
            continue
            
        protein1, pos1, protein2, pos2 = parts
        pos1, pos2 = int(pos1), int(pos2)
        
        match1 = find_matching_fragment(protein1, pos1, fragment_lookup)
        match2 = find_matching_fragment(protein2, pos2, fragment_lookup)
        
        if match1 and match2:
            real_pair = f"{match1[0]}:{match1[1]}:{match2[0]}:{match2[1]}"
            mapping[pair] = real_pair
    
    return mapping

def process_counts_files(input_dir: str, output_dir: str, real_fragments_list: List[str]):
    """Process all .counts files using optimized methods."""
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    # Build efficient fragment lookup
    print("Building fragment lookup...")
    fragment_lookup = build_fragment_lookup(real_fragments_list)
    
    # Process each .counts file
    for filename in os.listdir(input_dir):
        if filename.endswith('.counts'):
            input_path = os.path.join(input_dir, filename)
            output_path = os.path.join(output_dir, filename.replace('.counts', '.consolidated.counts'))
            
            print(f"Processing {filename}...")
            
            # Read file
            df = pd.read_csv(input_path, sep='\t')
            
            # Create fragment mapping once for all rows
            mapping = process_fragment_pairs(df, fragment_lookup)
            
            # Create result DataFrame efficiently
            if mapping:
                # Create a mask for rows that have matches
                mask = df['ID'].isin(mapping.keys())
                
                # Group and aggregate
                df['mapped_ID'] = df['ID'].map(mapping)
                result = df[mask].groupby('mapped_ID')[df.columns[1:]].sum()
                result.index.name = 'ID'
                
                # Save results
                result.to_csv(output_path, sep='\t')
                print(f"Saved consolidated results to {output_path}")
            else:
                print(f"No matching fragments found in {filename}")

def main():
    input_dir = "../CountsFiles/AllCounts"
    output_dir = "../CountsFiles/AllCounts/Consolidated"
    
    # Get real fragments
    real_fragments_list = fragments("sorted_genomeFrags.txt")
    
    # Process all counts files
    process_counts_files(input_dir, output_dir, real_fragments_list)

if __name__ == "__main__":
    main()