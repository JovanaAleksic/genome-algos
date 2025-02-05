import pandas as pd
import numpy as np
from collections import defaultdict

def parse_fragment(frag):
    """Parse fragment string into protein name and start position"""
    if ':' not in frag:
        return None, None
    protein, pos = frag.split(':')
    return protein, int(pos)

def calculate_overlap(start1, start2, fragment_length=140):
    """Calculate overlap between two fragments"""
    end1 = start1 + fragment_length
    end2 = start2 + fragment_length
    overlap_start = max(start1, start2)
    overlap_end = min(end1, end2)
    overlap = max(0, overlap_end - overlap_start)
    return overlap

def calculate_overlap_percentage(overlap, fragment_length=140):
    """Calculate overlap percentage only if there's at least 1 amino acid overlap"""
    if overlap >= 1:  # At least 1 amino acid overlap
        return (overlap / fragment_length) * 100
    return np.nan

def create_protein_pair_index(df):
    """Create an index of protein pairs to their fragment positions"""
    pair_index = defaultdict(list)
    
    # Pre-parse all fragments
    for idx, row in df.iterrows():
        protein1, start1 = parse_fragment(row['Fragment1'])
        protein2, start2 = parse_fragment(row['Fragment2'])
        if protein1 is not None and protein2 is not None:
            # Store both forward and reverse pairs
            pair_index[(protein1, protein2)].append((start1, start2, idx))
            if protein1 != protein2:  # Only add reverse if proteins are different
                pair_index[(protein2, protein1)].append((start2, start1, idx))
    
    return pair_index

def calculate_max_overlap_with_fragments(row_idx, protein1, start1, protein2, start2, pair_index):
    """Calculate maximum average overlap and return the matching fragments"""
    max_avg_overlap = 0
    best_frag1 = None
    best_frag2 = None
    
    # Check for same protein interactions
    if protein1 == protein2:
        overlap = calculate_overlap(start1, start2)
        overlap_perc = calculate_overlap_percentage(overlap)
        if not np.isnan(overlap_perc):  # Only consider valid overlaps
            max_avg_overlap = overlap_perc
            best_frag1 = f"{protein1}:{start1}"
            best_frag2 = f"{protein2}:{start2}"
    
    # Check other pairs with same proteins
    for other_start1, other_start2, other_idx in pair_index.get((protein1, protein2), []):
        if other_idx == row_idx:
            continue
        
        overlap1 = calculate_overlap(start1, other_start1)
        overlap2 = calculate_overlap(start2, other_start2)
        
        # Only calculate percentage if both fragments have at least 1 amino acid overlap
        if overlap1 >= 1 and overlap2 >= 1:
            overlap1_perc = calculate_overlap_percentage(overlap1)
            overlap2_perc = calculate_overlap_percentage(overlap2)
            avg_overlap = (overlap1_perc + overlap2_perc) / 2  # Take average of both overlaps
            
            if avg_overlap > max_avg_overlap:
                max_avg_overlap = avg_overlap
                best_frag1 = f"{protein1}:{other_start1}"
                best_frag2 = f"{protein2}:{other_start2}"
    
    return max_avg_overlap if max_avg_overlap > 0 else np.nan, best_frag1, best_frag2

def add_proximity_scores(df):
    """Add proximity scores and matching fragments to the dataframe"""
    # Create index of protein pairs
    pair_index = create_protein_pair_index(df)
    
    # Calculate scores and track matching fragments
    proximity_scores = []
    matching_frag1 = []
    matching_frag2 = []
    
    for idx, row in df.iterrows():
        protein1, start1 = parse_fragment(row['Fragment1'])
        protein2, start2 = parse_fragment(row['Fragment2'])
        
        if protein1 is None or protein2 is None:
            proximity_scores.append(np.nan)
            matching_frag1.append(None)
            matching_frag2.append(None)
            continue
        
        score, frag1, frag2 = calculate_max_overlap_with_fragments(
            idx, protein1, start1, protein2, start2, pair_index
        )
        
        proximity_scores.append(score)
        matching_frag1.append(frag1)
        matching_frag2.append(frag2)
    
    df['proximity_score'] = proximity_scores
    df['matching_fragment1'] = matching_frag1
    df['matching_fragment2'] = matching_frag2
    return df

def fragments(filename):
    fragments = []
    with open(filename, 'r') as f:
        next(f)  # Skip header line 'fragment_id'
        for line in f:
            fragments.append(line.strip())
    return fragments

# Main execution
real_fragments = fragments("sorted_genomeFrags.txt")
fs_fragments = fragments("../Results/fs_fragments_full.txt")
df = pd.read_csv('../Results/fragFrag_logFC1.5FDR0.05.txt', sep='\t')
df = df[df['Fragment1'].isin(real_fragments) & df['Fragment2'].isin(real_fragments)]
df = df[~(df['Fragment1'].isin(fs_fragments) & df['Fragment2'].isin(fs_fragments))]
df = add_proximity_scores(df)
print(df)
df.to_csv("../Results/fragFrag_logFC1.5FDR0.05_wproxperc.txt", index=False, sep="\t")