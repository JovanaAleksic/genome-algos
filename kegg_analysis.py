import requests
import pandas as pd
import numpy as np
from collections import defaultdict
from scipy import stats
import json
import time

def fetch_uniprot_from_refseq(protein_id):
    """Convert RefSeq ID to UniProt ID"""
    base_url = "https://rest.uniprot.org/uniprotkb/search?query="
    query = f"xref:RefSeq:{protein_id}"
    try:
        response = requests.get(f"{base_url}{query}&format=json")
        data = response.json()
        if 'results' in data and data['results']:
            return data['results'][0]['primaryAccession']
    except Exception as e:
        print(f"Error converting RefSeq ID {protein_id}: {str(e)}")
    return None

def get_kegg_data(uniprot_id):
    try:
        # First get KEGG ID from UniProt
        url = f"https://rest.kegg.jp/conv/hpy/{uniprot_id}"
        response = requests.get(url)
        time.sleep(1)  
        
        if response.status_code == 200:
            lines = response.text.strip().split('\n')
            if lines and lines[0]: 
                kegg_id = lines[0].split('\t')[1] if '\t' in lines[0] else None
                if kegg_id:
                    url = f"https://rest.kegg.jp/get/{kegg_id}"
                    response = requests.get(url)
                    time.sleep(1) 
                    
                    if response.status_code == 200:
                        pathways = []
                        for line in response.text.split('\n'):
                            if line.startswith('PATHWAY'):
                                pathway = line.split()[1]
                                pathways.append(pathway)
                        return pathways
    except Exception as e:
        print(f"Error fetching KEGG data for {uniprot_id}: {str(e)}")
    return []

def analyze_kegg_pathways(input_file):
    df = pd.read_csv(input_file, sep='\t')

    protein_pathways = {} 
    pathway_proteins = defaultdict(set)
    pathway_interactions = defaultdict(list)  
    all_proteins = set()  
    
    total_pairs = len(df)
    for idx, row in df.iterrows():
        print(f"Processing protein pair {idx + 1}/{total_pairs}")
        
        for protein_id in [row['Protein1'], row['Protein2']]:
            all_proteins.add(protein_id)
            
            if protein_id not in protein_pathways:
                # Convert RefSeq to UniProt
                uniprot_id = fetch_uniprot_from_refseq(protein_id)
                if uniprot_id:
                    # Get KEGG pathways
                    pathways = get_kegg_data(uniprot_id)
                    protein_pathways[protein_id] = pathways
                    
                    # Update pathway-protein mappings
                    for pathway in pathways:
                        pathway_proteins[pathway].add(protein_id)
    
        # Record interaction strengths for pathways
        protein1_pathways = protein_pathways.get(row['Protein1'], [])
        protein2_pathways = protein_pathways.get(row['Protein2'], [])
        
        for pathway in set(protein1_pathways + protein2_pathways):
            pathway_interactions[pathway].append({
                'strength': row['logFCmax'],
                'fdr': row['FDRmin'],
                'proteins': (row['Protein1'], row['Protein2'])
            })
    
    return perform_enrichment_analysis(pathway_proteins, pathway_interactions, all_proteins)

def perform_enrichment_analysis(pathway_proteins, pathway_interactions, all_proteins):
    """Calculate pathway enrichment statistics"""
    total_proteins = len(all_proteins)
    enrichment_results = []
    
    for pathway, proteins in pathway_proteins.items():
        pathway_size = len(proteins)
        num_interactions = len(pathway_interactions[pathway])
        
        # Calculate average interaction strength
        interaction_strengths = [inter['strength'] for inter in pathway_interactions[pathway]]
        avg_strength = np.mean(interaction_strengths) if interaction_strengths else 0
        
        # Hypergeometric test for enrichment
        M = total_proteins
        n = pathway_size
        N = len(all_proteins)
        k = num_interactions
        
        # Calculate p-value
        pvalue = stats.hypergeom.sf(k-1, M, n, N)
        
        enrichment_results.append({
            'pathway': pathway,
            'num_proteins': pathway_size,
            'num_interactions': num_interactions,
            'avg_interaction_strength': avg_strength,
            'pvalue': pvalue,
            'fdr': 1.0  # Initialize FDR to 1.0
        })
    
    # Sort by p-value and perform FDR correction
    if enrichment_results:
        enrichment_results.sort(key=lambda x: x['pvalue'])
        pvalues = [result['pvalue'] for result in enrichment_results]
        _, fdr_values, _, _ = stats.multipletests(pvalues, method='fdr_bh')
        for result, fdr in zip(enrichment_results, fdr_values):
            result['fdr'] = fdr
    
    return enrichment_results

def generate_pathway_reports(enrichment_results, output_prefix):
    """Generate analysis reports"""
    if not enrichment_results:
        print("No pathway enrichment results found.")
        with open(f"{output_prefix}_pathway_summary.txt", 'w') as f:
            f.write("KEGG Pathway Enrichment Analysis Summary\n")
            f.write("=====================================\n\n")
            f.write("No significant pathway enrichment found.\n")
        return
    
    results_df = pd.DataFrame(enrichment_results)
    results_df.to_csv(f"{output_prefix}_pathway_enrichment.csv", index=False)
    
    with open(f"{output_prefix}_pathway_summary.txt", 'w') as f:
        f.write("KEGG Pathway Enrichment Analysis Summary\n")
        f.write("=====================================\n\n")
        
        if len(results_df) > 0:
            f.write("Top 10 Enriched Pathways:\n")
            f.write("-----------------------\n")
            for _, row in results_df.head(10).iterrows():
                f.write(f"\nPathway: {row['pathway']}\n")
                f.write(f"Number of proteins: {row['num_proteins']}\n")
                f.write(f"Number of interactions: {row['num_interactions']}\n")
                f.write(f"Average interaction strength: {row['avg_interaction_strength']:.2f}\n")
                f.write(f"FDR-corrected p-value: {row['fdr']:.2e}\n")

            f.write("\nOverall Statistics:\n")
            f.write("-----------------\n")
            f.write(f"Total pathways analyzed: {len(results_df)}\n")
            significant_pathways = sum(results_df['fdr'] < 0.05)
            f.write(f"Pathways with FDR < 0.05: {significant_pathways}\n")
            f.write(f"Average proteins per pathway: {results_df['num_proteins'].mean():.1f}\n")
        else:
            f.write("No pathways found in the analysis.\n")

def main():
    input_file = "../Results/filtered_protProt_table.txt"
    output_prefix = "../Results/kegg_analysis"
    
    print("Starting KEGG pathway enrichment analysis...")
    enrichment_results = analyze_kegg_pathways(input_file)
    
    print("Generating reports...")
    generate_pathway_reports(enrichment_results, output_prefix)
    
    print("Analysis complete! Output files have been generated:")
    print(f"1. {output_prefix}_pathway_enrichment.csv")
    print(f"2. {output_prefix}_pathway_summary.txt")

if __name__ == "__main__":
    main()