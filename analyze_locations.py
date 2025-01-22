import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import networkx as nx
from collections import Counter
import numpy as np
from scipy.cluster import hierarchy
from scipy.spatial.distance import pdist
import itertools

def create_location_matrix(protein_locations, all_locations):
    """Create a binary matrix of proteins vs locations"""
    matrix = []
    proteins = []
    for protein, locs in protein_locations.items():
        row = [1 if loc in locs else 0 for loc in all_locations]
        matrix.append(row)
        proteins.append(protein)
    return pd.DataFrame(matrix, index=proteins, columns=all_locations)

def analyze_locations(file_path):
    # Read the TSV file
    df = pd.read_csv(file_path, sep='\t')
    
    # Filter for MDM2 interactions and valid locations
    mask = (df['GeneName2'] == 'MDM2') & (df['Location'].notna()) & \
           (df['Location'] != "Not MDM2 interaction") & \
           (df['Location'] != "No location data found")
    filtered_df = df[mask]
    
    # Process locations
    protein_locations = {}
    location_counts = Counter()
    all_locations_set = set()
    
    for _, row in filtered_df.iterrows():
        if pd.notna(row['Location']):
            locations = [loc.strip() for loc in row['Location'].split(';')]
            protein_locations[row['GeneName1']] = locations
            location_counts.update(locations)
            all_locations_set.update(locations)
    
    all_locations = sorted(list(all_locations_set))
    
    # Create the binary location matrix
    location_matrix = create_location_matrix(protein_locations, all_locations)
    
    # 1. Most common locations plot
    plt.figure(figsize=(15, 8))
    most_common = dict(location_counts.most_common(15))
    sns.barplot(x=list(most_common.values()), y=list(most_common.keys()))
    plt.title('Top 15 Most Common Cellular Locations Among MDM2 Interactors')
    plt.xlabel('Count')
    plt.ylabel('Location')
    plt.tight_layout()
    plt.show()
    
    # 2. Location correlation analysis
    plt.figure(figsize=(15, 12))
    # Get top 30 locations for correlation analysis to keep visualization manageable
    top_locations = [loc for loc, _ in location_counts.most_common(30)]
    correlation_matrix = location_matrix[top_locations].corr()
    
    sns.heatmap(correlation_matrix, cmap='coolwarm', center=0,
                xticklabels=True, yticklabels=True)
    plt.title('Location Correlation Matrix (Top 30 Locations)')
    plt.xticks(rotation=90)
    plt.yticks(rotation=0)
    plt.tight_layout()
    plt.show()
    
    # 3. Hierarchical clustering of proteins
    # Limit to proteins with at least 2 locations for meaningful clustering
    proteins_multi_loc = {p: locs for p, locs in protein_locations.items() if len(locs) >= 2}
    if len(proteins_multi_loc) > 3:  # Ensure we have enough data to cluster
        matrix_multi = create_location_matrix(proteins_multi_loc, all_locations)
        
        plt.figure(figsize=(15, 10))
        Z = hierarchy.linkage(pdist(matrix_multi), method='ward')
        hierarchy.dendrogram(Z, labels=matrix_multi.index, leaf_rotation=90)
        plt.title('Hierarchical Clustering of Proteins Based on Location Patterns\n(Proteins with ≥2 locations)')
        plt.xlabel('Proteins')
        plt.ylabel('Distance')
        plt.tight_layout()
        plt.show()


    
    # Print summary statistics
    print("\nSummary Statistics:")
    print(f"Total proteins analyzed: {len(protein_locations)}")
    print(f"Total unique locations: {len(location_counts)}")
    print(f"Average locations per protein: {np.mean([len(locs) for locs in protein_locations.values()]):.2f}")
    print(f"\nSubcellular Distribution:")
    print(f"Nuclear proteins: {nuclear_count}")
    print(f"Cytoplasmic proteins: {cytoplasmic_count}")
    print(f"Proteins with both nuclear and cytoplasmic localization: {both_count}")
    
    # Print top 10 most common locations
    print("\nTop 10 most common locations:")
    for loc, count in location_counts.most_common(10):
        print(f"{loc}: {count}")
    
    return protein_locations, location_counts, location_matrix

if __name__ == "__main__":
    file_path = "mdm2_APIDinteractions_with_locations.tsv"
    protein_locations, location_counts, location_matrix = analyze_locations(file_path)