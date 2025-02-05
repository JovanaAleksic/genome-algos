import networkx as nx
import numpy as np
import pandas as pd
from collections import Counter

def load_network(filename):
    # Read your tab-separated interaction data
    df = pd.read_csv(filename, sep='\t')
    G = nx.Graph()
    for _, row in df.iterrows():
        G.add_edge(row['Protein1'], row['Protein2'])
    return G

def get_network_properties(G):
    props = {
        'clustering': nx.average_clustering(G),
        'density': nx.density(G),
        'avg_degree': np.mean([d for n, d in G.degree()])
    }
    # Add average path length only if network is connected
    if nx.is_connected(G):
        props['avg_path_length'] = nx.average_shortest_path_length(G)
    return props

def generate_random_networks(G, n_random=1000):
    degree_sequence = [d for n, d in G.degree()]
    random_properties = []
    
    for i in range(n_random):
        G_random = nx.configuration_model(degree_sequence)
        G_random = nx.Graph(G_random)  # Remove parallel edges and self-loops
        random_properties.append(get_network_properties(G_random))
    
    return random_properties

def calculate_zscores(orig_props, random_props):
    z_scores = {}
    for metric in orig_props.keys():
        random_values = [r[metric] for r in random_props if metric in r]
        if not random_values:
            continue
        mean = np.mean(random_values)
        std = np.std(random_values)
        if std > 0:
            z_scores[metric] = (orig_props[metric] - mean) / std
    return z_scores

# Main analysis
G = load_network('../Results/protProt_filtFragPairs2Libs_logFC1.5FDR0.05_woutSticky.txt')
orig_properties = get_network_properties(G)
random_properties = generate_random_networks(G, n_random=1000)
z_scores = calculate_zscores(orig_properties, random_properties)

print("Original network properties:")
for k, v in orig_properties.items():
    print(f"{k}: {v:.4f}")

print("\nZ-scores (compared to random networks):")
for k, v in z_scores.items():
    print(f"{k}: {v:.4f}")  # Changed 'z' to 'v'