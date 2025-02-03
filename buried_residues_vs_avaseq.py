import pandas as pd
import numpy as np

def parse_bsa_report(bsa_file):
    """Parse the buried surface area report content."""
    buried_residues = {}
    current_complex = None
    current_rank = None
    
    with open(bsa_file, 'r') as f:
        for line in f:
            if line.startswith('Complex:'):
                current_complex = line.split(':')[1].strip()
                buried_residues[current_complex] = {}
                
            elif line.startswith('Rank'):
                rank_num = int(line.split(':')[0].split()[1])
                current_rank = rank_num
                buried_residues[current_complex][rank_num] = {'A': set(), 'B': set(), 'C': set()}
                
            elif 'Chain' in line and ':' in line:
                parts = line.split(':')
                chain = parts[0].strip().split()[-1]
                residues = parts[1].strip()
                # Parse residue numbers
                res_nums = [int(x) for x in residues.split(',')]
                buried_residues[current_complex][current_rank][chain].update(res_nums)
                
    return buried_residues

def analyze_fragment_overlap(fragment_start, is_mdm2, buried_residues):
    """Calculate what fraction of the fragment residues are buried."""
    if is_mdm2:
        # For MDM2, fragment extends from start to residue 491
        fragment = set(range(fragment_start, 492))  # Include 491
        fragment_length = 492 - fragment_start
    else:
        # For other proteins, use 120 residue fragments
        fragment = set(range(fragment_start, fragment_start + 120))
        fragment_length = 120

    if not buried_residues:
        return 0
    
    return len(fragment.intersection(buried_residues)) / fragment_length

def get_chain_mapping(buried_data_complex):
    """
    Determine the correct chain mapping for a complex.
    MDM2 is always the second chain (B or C) in the buried surface predictions.
    Returns tuple of (partner_chain, mdm2_chain)
    """
    # Check if complex has A/B or B/C chains
    if 'B' in buried_data_complex[0] and 'A' in buried_data_complex[0]:  # A/B combination
        mdm2_chain = 'B'
        partner_chain = 'A'
    else:  # B/C combination
        mdm2_chain = 'C'
        partner_chain = 'B'
    
    return (partner_chain, mdm2_chain)

def main():
    # Read files
    df = pd.read_csv('../Results/ffi_libs2+.txt')
    
    # Parse BSA report
    buried_data = parse_bsa_report('bsa_analysis_report.txt')
    
    # Analyze each interaction
    match_percentages = []
    best_ranks = []
    
    for _, row in df.iterrows():
        # Get proteins and fragments
        protein1, protein2 = row['protein1'], row['protein2']
        frag1_start = int(row['frag1'])
        frag2_start = int(row['frag2'])
        
        # Determine which protein is MDM2 and create proper complex name
        mdm2_is_first = (protein1 == 'MDM2')
        if mdm2_is_first:
            # Swap order for complex name since MDM2 should be second
            complex_name = f"{protein2}_{protein1}"
        else:
            complex_name = f"{protein1}_{protein2}"
        
        best_match_score = 0
        best_match_rank = None
        
        # Check if complex exists in buried data
        if complex_name in buried_data:
            # Get the correct chain mapping
            partner_chain, mdm2_chain = get_chain_mapping(buried_data[complex_name])
            
            # Check each rank prediction
            for rank, chains in buried_data[complex_name].items():
                # If MDM2 is first protein, swap the fragment analysis
                if mdm2_is_first:
                    mdm2_score = analyze_fragment_overlap(frag1_start, True, chains[mdm2_chain])
                    partner_score = analyze_fragment_overlap(frag2_start, False, chains[partner_chain])
                else:
                    partner_score = analyze_fragment_overlap(frag1_start, False, chains[partner_chain])
                    mdm2_score = analyze_fragment_overlap(frag2_start, True, chains[mdm2_chain])
                
                avg_score = (partner_score + mdm2_score) / 2
                
                if avg_score > best_match_score:
                    best_match_score = avg_score
                    best_match_rank = rank
        
        # Store results, leave blank for zero percentages
        percentage = round(best_match_score * 100, 2)
        match_percentages.append(percentage if percentage > 0 else None)
        best_ranks.append(best_match_rank)
    
    # Add new columns to dataframe
    df['match_percentage'] = match_percentages
    
    # Convert best_match_rank to integer type (while preserving NaN)
    df['best_match_rank'] = pd.Series(best_ranks, dtype='Int64')
    
    # Save results to CSV
    df.to_csv('../Results/ffi_wburied.csv', index=False, na_rep='', float_format='%.2f')
    print("Analysis complete. First few rows of results:")
    print(df.head())
    return df

if __name__ == "__main__":
    result_df = main()