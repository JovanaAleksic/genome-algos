import pandas as pd
import numpy as np
import os
from pathlib import Path
from itertools import combinations

def transform_fragment_name(name):

    parts = name.split("_")
    position = int(int(parts[2])/3) + 1
    return f"{parts[0]}_{parts[1]}:{position}"

def normalize_pair(pair):

    # Split the pair into two fragments
    parts = pair.split(':')

    protein1, start1, protein2, start2 = parts
    frag1 = protein1 + ":" + start1
    frag2 = protein2 + ":" + start2
    
    # Sort fragments alphabetically to get consistent orientation
    if frag1 < frag2:
        return f"{frag1}:{frag2}"
    else:
        return f"{frag2}:{frag1}"

def generate_all_possible_pairs(fragments_file):

    # Read fragments file
    fragments_df = pd.read_csv(fragments_file)
    
    # Filter for pools 1-5
    valid_pools = list(range(1, 6))  # pools 1-5
    fragments_df['Pool'] = fragments_df['Pool'].astype(int)
    fragments_df = fragments_df[fragments_df['Pool'].isin(valid_pools)]
    
    # Transform fragment names
    transformed_names = fragments_df['Name'].apply(transform_fragment_name).unique()
    
    # Generate all possible pairs (only one orientation per pair)
    all_pairs = set()
    for pair in combinations(sorted(transformed_names), 2):
        # Add only normalized version of the pair
        all_pairs.add(':'.join(sorted(pair)))
    
    print(f"Generated {len(all_pairs)} possible pairs from {len(transformed_names)} fragments")
    return all_pairs

def analyze_fragment_pairs_across_files(counts_directory, all_possible_pairs):

    # Dictionary to store how many files each pair appears in
    pair_appearances = {}
    invalid_pairs = set()
    
    # Get list of all .counts files
    counts_files = list(Path(counts_directory).glob("*.counts"))
    total_files = len(counts_files)
    
    print(f"Found {total_files} count files to analyze...")
    
    # Process each counts file
    for counts_file in counts_files:
        try:
            print(f"Processing {counts_file.name}...")
            
            # Read the file
            df = pd.read_csv(counts_file, sep='\t')
            
            # Get non-zero counts
            counts_matrix = df.iloc[:, 1:] # exclude ID column
            non_zero_rows = counts_matrix.any(axis=1)
            present_pairs = df.loc[non_zero_rows, 'ID'].tolist()
            
            # Normalize and update appearance count for each present pair
            for pair in present_pairs:
                try:
                    normalized_pair = normalize_pair(pair)
                    # Only count the pair if it's in our valid combinations
                    if normalized_pair in all_possible_pairs:
                        pair_appearances[normalized_pair] = pair_appearances.get(normalized_pair, 0) + 1
                    else:
                        invalid_pairs.add(normalized_pair)
                except Exception as e:
                    print(f"Warning: Could not process pair {pair} in file {counts_file.name}")
                    continue
                
        except Exception as e:
            print(f"Error processing file {counts_file}: {str(e)}")
            continue
    
    # Analyze results
    pairs_in_multiple = {pair: count for pair, count in pair_appearances.items() 
                        if count >= 2}
    
    results = {
        'total_files_processed': total_files,
        'total_possible_pairs': len(all_possible_pairs),
        'total_observed_valid_pairs': len(pair_appearances),
        'total_invalid_pairs': len(invalid_pairs),
        'pairs_in_multiple_files': len(pairs_in_multiple),
        'coverage_percentage': (len(pairs_in_multiple) / len(all_possible_pairs)) * 100 if all_possible_pairs else 0,
        'pair_counts': pairs_in_multiple,
        'invalid_pairs': invalid_pairs
    }
    
    return results

def print_results(results):

    print("\nAnalysis Results:")
    print("-" * 50)
    print(f"Total files processed: {results['total_files_processed']}")
    print(f"Total possible pairs (from pools 1-5): {results['total_possible_pairs']}")
    print(f"Total observed valid pairs: {results['total_observed_valid_pairs']}")
    print(f"Total invalid pairs found: {results['total_invalid_pairs']}")
    print(f"Fragment pairs appearing in multiple files: {results['pairs_in_multiple_files']}")
    print(f"Coverage percentage (pairs in ≥2 files): {results['coverage_percentage']:.2f}%")
    
    print("\nTop 10 most frequent valid pairs:")
    print("-" * 50)
    sorted_pairs = sorted(results['pair_counts'].items(), 
                         key=lambda x: x[1], 
                         reverse=True)
    
    for pair, count in sorted_pairs[:10]:
        print(f"Pair: {pair}")
        print(f"Appears in {count} files ({(count/results['total_files_processed'])*100:.1f}% of files)")
        print("-" * 30)
    
    print("\nFirst 10 invalid pairs found (if any):")
    print("-" * 50)
    for pair in list(results['invalid_pairs'])[:10]:
        print(f"Invalid pair: {pair}")

def save_results(results, output_file="fragment_pairs_analysis.txt"):

    with open(output_file, 'w') as f:
        f.write("Fragment Pairs Analysis Results (Pools 1-5)\n")
        f.write("=" * 50 + "\n\n")
        
        
        f.write("Summary Statistics:\n")
        f.write("----------------------------------------------" + "\n")
        f.write(f"Total files processed: {results['total_files_processed']}\n")
        f.write(f"Total possible pairs: {results['total_possible_pairs']}\n")
        f.write(f"Total observed valid pairs: {results['total_observed_valid_pairs']}\n")
        f.write(f"Total invalid pairs found: {results['total_invalid_pairs']}\n")
        f.write(f"Pairs in multiple files: {results['pairs_in_multiple_files']}\n")
        f.write(f"Coverage percentage: {results['coverage_percentage']:.2f}%\n\n")
        
        f.write("Detailed Valid Pair Counts:\n")
        f.write("-" * 30 + "\n")
        sorted_pairs = sorted(results['pair_counts'].items(), 
                            key=lambda x: x[1], 
                            reverse=True)
        
        for pair, count in sorted_pairs:
            f.write(f"{pair}\t{count}\t{(count/results['total_files_processed'])*100:.1f}%\n")
            
        f.write("\nInvalid Pairs Found:\n")
        f.write("-" * 30 + "\n")
        for pair in sorted(results['invalid_pairs']):
            f.write(f"{pair}\n")

def main():

    fragments_file = "../FragmentList.csv"
    counts_directory = "../CountsFiles/AllCounts/"
    
    print("Generating all possible fragment pairs")
    all_possible_pairs = generate_all_possible_pairs(fragments_file)
    
    results = analyze_fragment_pairs_across_files(counts_directory, all_possible_pairs)
    
    print_results(results)
    save_results(results)

if __name__ == "__main__":
    main()