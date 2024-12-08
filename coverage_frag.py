import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict
import os
from concurrent.futures import ProcessPoolExecutor
import pandas as pd

def parse_fragments_list(file_obj):
    return [line.strip() for line in file_obj if line.strip()]

def process_counts_file(args):
    filename, fragments_list = args
    local_pairings = defaultdict(set)
    
    # Read file
    try:
        print(filename)
        # Read only the necessary columns (first column and count columns)
        df = pd.read_csv(filename, sep='\t')
        
        # Filter rows where sum of counts >= 10
        df = df[df.iloc[:, 1:].sum(axis=1) >= 10]
        
        # Process filtered rows
        for _, row in df.iterrows():
            frag_pair = row.iloc[0].split(':')
            if len(frag_pair) < 3:
                continue
                
            frag1 = f"{frag_pair[0]}:{frag_pair[1]}"
            frag2 = f"{frag_pair[2]}:{frag_pair[3]}"
            
            if frag1 in fragments_list and frag2 in fragments_list:
                key = (min(frag1, frag2), max(frag1, frag2))
                local_pairings[key].add(os.path.basename(filename))
        
        return dict(local_pairings)

    except Exception as e:
        print(f"Error processing {filename}: {str(e)}")
        return {}

def analyze_fragment_pairings(fragments_list, counts_dir):
    # fragments_set = set(fragments_list)
    
    # Get list of count files
    count_files = [os.path.join(counts_dir, f) for f in os.listdir(counts_dir) 
                  if f.endswith('.counts')]
    
    # Process files in parallel
    args = [(f, fragments_list) for f in count_files]
    all_pairings = defaultdict(set)
    
    with ProcessPoolExecutor() as executor:
        results = executor.map(process_counts_file, args)
        
        # Combine results
        for file_pairings in results:
            for pair, files in file_pairings.items():
                all_pairings[pair].update(files)
    
    # Count significant pairings for each fragment
    fragment_counts = defaultdict(int)
    for (frag1, frag2), files in all_pairings.items():
        if len(files) >= 1:
            fragment_counts[frag1] += 1
            fragment_counts[frag2] += 1
    
    return fragment_counts

def plot_fragment_pairings(fragments_list, fragment_counts):
    plt.figure(figsize=(20, 8))
    
    x = np.arange(len(fragments_list))
    y = [fragment_counts.get(frag, 0) for frag in fragments_list]
    
    plt.bar(x, y, width=1.0)
    
    plt.title('Number of Significant Fragment Pairings (≥2 libs, counts sum ≥10)', fontsize=14, pad=20)
    plt.xlabel('Fragments (in genome order)', fontsize=12)
    plt.ylabel('Number of Fragment Pairings', fontsize=12)

    zero_pairing_indices = [i for i, count in enumerate(y) if count == 0]
    zero_pairing_labels = [fragments_list[i] for i in zero_pairing_indices]
    
    plt.xticks(zero_pairing_indices, zero_pairing_labels, rotation=45, ha='right')
    plt.grid(True, axis='y', linestyle='--', alpha=0.7)
    
    plt.tight_layout()
    plt.savefig('fragment_pairings_plotL.png', dpi=300, bbox_inches='tight')
    
    print(f"Maximum pairings: {max(y)}")
    print(f"Average pairings: {np.mean(y):.2f}")
    print(f"Number of fragments with at least one pairing: {sum(1 for count in y if count > 0)}")

def main():
    print("Reading fragments list...")
    with open('../sorted_genomeFrags.txt', 'r') as f:
        fragments_list = parse_fragments_list(f)
    
    fragment_counts = analyze_fragment_pairings(fragments_list, '../../CountsFiles/AllCounts/')
    
    plot_fragment_pairings(list(fragments_list), fragment_counts)
    plt.show()

if __name__ == "__main__":
    main()