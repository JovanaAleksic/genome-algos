import requests
import pandas as pd
import statistics

# def get_string_interactions(protein, species=9606):
#     """
#     Get protein interactions from STRING database
#     """
#     string_api_url = "https://string-db.org/api"
#     output_format = "json"
#     method = "network"
    
#     request_url = "/".join([string_api_url, output_format, method])
    
#     params = {
#         "identifiers": protein,
#         "species": species,
#         "required_score": 400,  # Medium confidence (0.4)
#         "network_type": "physical"
#     }
    
#     try:
#         response = requests.post(request_url, data=params)
#         return response.json()
#     except:
#         return None

# def analyze_interactions(protein_list):
#     """
#     Analyze interactions for a list of proteins
#     """
#     # Store interaction counts
#     interaction_counts = {}
#     failed_proteins = []
    
#     # Get interactions for each protein
#     for protein in protein_list:
#         interactions = get_string_interactions(protein)
#         if interactions:
#             interaction_counts[protein] = len(interactions)
#         else:
#             failed_proteins.append(protein)
            
#     # Calculate statistics
#     if interaction_counts:
#         avg_interactions = statistics.mean(interaction_counts.values())
#         median_interactions = statistics.median(interaction_counts.values())
        
#         # Create DataFrame for results
#         results_df = pd.DataFrame.from_dict(interaction_counts, orient='index', 
#                                           columns=['interaction_count'])
#         results_df.sort_values('interaction_count', ascending=False, inplace=True)
        
#         return {
#             'average_interactions': avg_interactions,
#             'median_interactions': median_interactions,
#             'detailed_results': results_df,
#             'failed_proteins': failed_proteins
#         }
#     return None

# # List of proteins to analyze
# proteins = [
#     "MCM10", "BCL2L1", "IFIT1", "C3orf38", "RCC1", "PSMD4", "GRAP2", "ORC2L", 
#     "CRADD", "PEX19", "LCP2", "EIF3E", "SLC22A15", "AKT1", "IGFBP4", "LMNB1", 
#     "SMAD4", "PPP6C", "PSMD5", "PDPK1", "IGF2", "LMNA", "SMAD1", "ZNF350", 
#     "PEX14", "SYCE1", "GTF2F2", "GADD45A", "DCP1A", "CASP2", "BDNF", "ARF1", 
#     "BAK1", "GTF2F1", "PPP3CA", "XIAP", "JUNB", "NTF4", "AFRIP2", "FOS", 
#     "FEN1", "PPP3R1", "CASP9", "BATF", "PEX3", "ERBB3", "GCDH", "CEBPG", 
#     "PCNA", "RAC1", "CASP7", "BAD", "ORC4L", "NRG1", "ZCCHC9", "RAN", "RAD23A"
# ]

# # Run the analysis
# results = analyze_interactions(proteins)

# # Print results
# if results:
#     print(f"\nAverage interactions for tested proteins: {results['average_interactions']:.2f}")
#     print(f"Median interactions for tested proteins: {results['median_interactions']:.2f}")
#     print("\nDetailed interaction counts:")
#     print(results['detailed_results'])
    
#     if results['failed_proteins']:
#         print("\nFailed to fetch data for these proteins:")
#         print(results['failed_proteins'])



def analyze_apid_interactions(file_path, proteins_of_interest):
    # Read the APID file with tab separation and custom handling
    try:
        # First try with tab separator
        df = pd.read_csv(file_path, sep='\t')
    except:
        try:
            # If tab doesn't work, try reading with whitespace separation
            df = pd.read_csv(file_path, delim_whitespace=True)
        except:
            # If both fail, try reading the file manually
            with open(file_path, 'r') as f:
                lines = f.readlines()
            
            # Parse header
            header = lines[0].strip().split()
            
            # Parse data lines
            data = []
            for line in lines[1:]:
                # Split on whitespace while preserving quoted strings
                parts = line.strip().split()
                if len(parts) >= 7:  # Ensure we have at least the essential columns
                    data.append(parts)
            
            # Create DataFrame
            df = pd.DataFrame(data, columns=header)
    
    # Convert proteins of interest to uppercase for consistent matching
    proteins_of_interest = [p.upper() for p in proteins_of_interest]
    
    # Count interactions for each protein (both as A and B)
    protein_counts_a = df['GeneName_A'].value_counts()
    protein_counts_b = df['GeneName_B'].value_counts()
    
    # Combine counts (add counts where protein appears as both A and B)
    total_counts = protein_counts_a.add(protein_counts_b, fill_value=0)
    
    # Calculate average interactions for all proteins
    overall_average = total_counts.mean()
    
    # Get counts for proteins of interest
    poi_counts = total_counts[total_counts.index.isin(proteins_of_interest)]
    poi_average = poi_counts.mean() if not poi_counts.empty else 0
    
    # Create detailed report for proteins of interest
    poi_details = pd.DataFrame({
        'Protein': poi_counts.index,
        'Interaction_Count': poi_counts.values
    }).sort_values('Interaction_Count', ascending=False)
    
    # Find proteins of interest that weren't found in the database
    missing_proteins = set(proteins_of_interest) - set(total_counts.index)
    
    # Calculate some additional statistics
    stats = {
        'overall_average': overall_average,
        'overall_median': total_counts.median(),
        'poi_average': poi_average,
        'poi_median': poi_counts.median() if not poi_counts.empty else 0,
        'poi_details': poi_details,
        'missing_proteins': missing_proteins,
        'total_proteins': len(total_counts),
        'total_interactions': len(df)
    }
    
    return stats

# List of proteins you tested
proteins_of_interest = [
    "MCM10", "BCL2L1", "IFIT1", "C3orf38", "RCC1", "PSMD4", "GRAP2", "ORC2", 
    "CRADD", "PEX19", "LCP2", "EIF3E", "SLC22A15", "AKT1", "IGFBP4", "LMNB1", 
    "SMAD4", "PPP6C", "PSMD5", "PDPK1", "IGF2", "LMNA", "SMAD1", "ZNF350", 
    "PEX14", "SYCE1", "GTF2F2", "GADD45A", "DCP1A", "CASP2", "BDNF", "ARF1", 
    "BAK1", "GTF2F1", "PPP3CA", "XIAP", "JUNB", "NTF4", "ARFIP2", "FOS", 
    "FEN1", "PPP3R1", "CASP9", "BATF", "PEX3", "ERBB3", "GCDH", "CEBPG", 
    "PCNA", "RAC1", "CASP7", "BAD", "ORC4", "NRG1", "ZCCHC9", "RAN", "RAD23A"
]

# Analyze the data
try:
    results = analyze_apid_interactions('../apid_human.txt', proteins_of_interest)
    
    # Print results
    print("\nDatabase Summary:")
    print(f"Total proteins in database: {results['total_proteins']:,}")
    print(f"Total interactions in database: {results['total_interactions']:,}")
    print(f"\nOverall Database Statistics:")
    print(f"Average interactions per protein: {results['overall_average']:.2f}")
    print(f"Median interactions per protein: {results['overall_median']:.2f}")
    print(f"\nTested Proteins Statistics:")
    print(f"Average interactions: {results['poi_average']:.2f}")
    print(f"Median interactions: {results['poi_median']:.2f}")
    
    print("\nDetailed interaction counts for tested proteins:")
    print(results['poi_details'])
    
    if results['missing_proteins']:
        print("\nWarning: The following proteins were not found in the database:")
        print(sorted(results['missing_proteins']))

except Exception as e:
    print(f"An error occurred: {str(e)}")
    print("Please check the file path and format.")
print(results['poi_details'])

if results['missing_proteins']:
    print("\nWarning: The following proteins were not found in the database:")
    print(sorted(results['missing_proteins']))