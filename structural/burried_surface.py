from Bio.PDB import *
import numpy as np
from typing import Dict, List, Tuple
import os
import pandas as pd
from collections import defaultdict

# Import all the existing functions from the original code
def calculate_residue_sasa(structure: Structure, chain_id: str = None) -> Dict[str, float]:
    """
    Calculate SASA for each residue in the structure or specific chain
    """
    sr = ShrakeRupley()
    sr.compute(structure, level="R")
    
    residue_sasa = {}
    for model in structure:
        for chain in model:
            if chain_id is None or chain.id == chain_id:
                for residue in chain:
                    res_id = f"{chain.id}_{residue.id[1]}_{residue.resname}"
                    residue_sasa[res_id] = residue.sasa
                    
    return residue_sasa

def identify_buried_residues(
    individual_sasa: Dict[str, float], 
    complex_sasa: Dict[str, float], 
    threshold: float = 1.0
) -> List[str]:
    """
    Identify buried residues based on SASA difference
    """
    buried_residues = []
    
    for res_id in individual_sasa:
        if res_id in complex_sasa:
            sasa_diff = individual_sasa[res_id] - complex_sasa[res_id]
            if sasa_diff > threshold:
                buried_residues.append(res_id)
                
    return buried_residues

def calculate_bsa(pdb_file: str) -> Tuple[float, Dict[str, List[str]]]:
    """
    Calculate buried surface area and identify buried residues in a protein complex
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('complex', pdb_file)
    
    chains = list(structure.get_chains())
    if len(chains) != 2:
        raise ValueError("Expected exactly 2 chains in the complex")
    
    sr = ShrakeRupley()
    sr.compute(structure, level="S")
    complex_area = structure.sasa
    
    complex_residue_sasa = calculate_residue_sasa(structure)
    
    individual_areas = []
    buried_residues = {}
    
    for chain in chains:
        tmp_structure = Structure.Structure('temp')
        tmp_model = Model.Model(0)
        tmp_structure.add(tmp_model)
        tmp_chain = Chain.Chain(chain.id)
        tmp_model.add(tmp_chain)
        
        for residue in chain:
            tmp_chain.add(residue.copy())
            
        sr.compute(tmp_structure, level="S")
        individual_areas.append(tmp_structure.sasa)
        
        individual_residue_sasa = calculate_residue_sasa(tmp_structure, chain.id)
        
        buried = identify_buried_residues(individual_residue_sasa, complex_residue_sasa)
        buried_residues[chain.id] = buried
    
    total_individual_area = sum(individual_areas)
    buried_surface_area = total_individual_area - complex_area
    
    return buried_surface_area, buried_residues

def get_residue_number(res_id: str) -> int:
    """Extract residue number from identifier"""
    chain, number, name = res_id.split('_')
    return int(number)

def process_complex(pdb_file: str) -> Dict:
    """
    Process a single complex and return its analysis results
    
    Args:
        pdb_file: Path to the PDB file
        
    Returns:
        Dictionary containing analysis results
    """
    try:
        bsa, buried_residues = calculate_bsa(pdb_file)
        
        # Get complex name and rank from file path
        dir_name = os.path.dirname(pdb_file)
        complex_name = os.path.basename(dir_name)
        rank = int(os.path.basename(pdb_file).split('_')[1].split('.')[0])
        
        result = {
            'complex_name': complex_name,
            'rank': rank,
            'bsa': bsa,
            'buried_residues': {}
        }
        
        for chain_id, residues in buried_residues.items():
            residue_numbers = [str(get_residue_number(r)) for r in sorted(residues)]
            result['buried_residues'][chain_id] = ','.join(residue_numbers)
            
        return result
        
    except Exception as e:
        print(f"Error processing {pdb_file}: {str(e)}")
        return None

def analyze_all_complexes(base_dir: str) -> pd.DataFrame:
    """
    Analyze all complexes in the given directory
    
    Args:
        base_dir: Base directory containing complex subdirectories
        
    Returns:
        DataFrame containing analysis results
    """
    results = []
    
    # Walk through all subdirectories
    for root, dirs, files in os.walk(base_dir):
        # Filter for ranked_*.pdb files
        pdb_files = [f for f in files if f.startswith('ranked_0') and f.endswith('.pdb')]
        
        for pdb_file in pdb_files:
            full_path = os.path.join(root, pdb_file)
            print("Working on: ", full_path)
            result = process_complex(full_path)
            
            if result:
                results.append(result)
    
    # Convert results to DataFrame
    df = pd.DataFrame(results)
    
    # Sort by complex name and rank
    if not df.empty:
        df = df.sort_values(['complex_name', 'rank'])
    
    return df

def generate_report(df: pd.DataFrame, output_file: str = 'bsa_analysis_report.txt'):
    """
    Generate a detailed report from the analysis results
    
    Args:
        df: DataFrame containing analysis results
        output_file: Path to output file
    """
    with open(output_file, 'w') as f:
        f.write("Buried Surface Area Analysis Report\n")
        f.write("=================================\n\n")
        
        # Group by complex
        for complex_name in df['complex_name'].unique():
            complex_data = df[df['complex_name'] == complex_name]
            
            f.write(f"Complex: {complex_name}\n")
            f.write("-" * (len(complex_name) + 9) + "\n")
            
            for _, row in complex_data.iterrows():
                f.write(f"\nRank {row['rank']}:\n")
                f.write(f"  Buried Surface Area: {row['bsa']:.2f} Å²\n")
                f.write("  Buried Residues:\n")
                
                buried_residues = eval(row['buried_residues']) if isinstance(row['buried_residues'], str) else row['buried_residues']
                for chain, residues in buried_residues.items():
                    f.write(f"    Chain {chain}: {residues}\n")
            
            f.write("\n" + "=" * 50 + "\n\n")
        
        # Add summary statistics
        f.write("\nSummary Statistics\n")
        f.write("=================\n")
        f.write(f"Total complexes analyzed: {len(df['complex_name'].unique())}\n")
        f.write(f"Total conformations analyzed: {len(df)}\n")
        f.write(f"Average BSA across all conformations: {df['bsa'].mean():.2f} Å²\n")
        f.write(f"Maximum BSA: {df['bsa'].max():.2f} Å²\n")
        f.write(f"Minimum BSA: {df['bsa'].min():.2f} Å²\n")

def main():
    """
    Main function to process all complexes and generate report
    """
    base_dir = "../Predictions/AF/"
    
    print("Analyzing complexes...")
    results_df = analyze_all_complexes(base_dir)
    
    if results_df.empty:
        print("No complexes found or all analyses failed.")
        return
    
    print("\nGenerating report...")
    generate_report(results_df)
    
    # Also save results as CSV for further analysis
    results_df.to_csv('bsa_analysis_results.csv', index=False)
    
    print("\nAnalysis complete!")
    print("Report saved as 'bsa_analysis_report.txt'")
    print("Detailed results saved as 'bsa_analysis_results.csv'")

if __name__ == "__main__":
    main()