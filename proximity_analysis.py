# import pandas as pd
# import numpy as np
# from scipy import stats
# from typing import List, Tuple, Dict
# import matplotlib.pyplot as plt

# class ProteinProximityAnalyzer:
#     def __init__(self, ppi_file: str, genome_order_file: str):
#         """
#         Initialize the analyzer with protein interaction and genome order data
        
#         Parameters:
#         -----------
#         ppi_file : str
#             Path to Excel file containing protein-protein interactions
#             Expected columns: Protein1, Protein2, plus additional metadata
#         genome_order_file : str
#             Path to file containing protein order in genome
#             Expected format: one column with protein IDs in genomic order
#         """
#         self.ppis = self._load_interactions(ppi_file)
#         self.protein_positions = self._load_genome_order(genome_order_file)
        
#     def _load_interactions(self, ppi_file: str) -> List[Tuple[str, str]]:
#         """Load protein-protein interactions from Excel file"""
#         # Read Excel file
#         df = pd.read_excel(ppi_file)
        
#         # Extract just the protein pairs
#         protein_pairs = list(zip(df['Protein1'], df['Protein2']))
        
#         # Store the full dataframe for later analysis if needed
#         self.full_ppi_data = df
        
#         return protein_pairs
    
#     def _load_genome_order(self, genome_file: str) -> Dict[str, int]:
#         """Load genomic positions and create position lookup dictionary"""
#         df = pd.read_csv(genome_file)
#         # Assuming single column with protein IDs
#         proteins = df.iloc[:, 0].tolist()
#         return {protein: pos for pos, protein in enumerate(proteins)}
    
#     def calculate_genomic_distances(self) -> Dict[str, List]:
#         """Calculate genomic distances between all interacting proteins"""
#         distances = []
#         pairs = []
#         logFCmax_values = []
        
#         for idx, row in self.full_ppi_data.iterrows():
#             prot1, prot2 = row['Protein1'], row['Protein2']
            
#             # Skip if either protein is not in the genome order list
#             if prot1 not in self.protein_positions or prot2 not in self.protein_positions:
#                 continue
                
#             pos1 = self.protein_positions[prot1]
#             pos2 = self.protein_positions[prot2]
            
#             # Calculate shortest distance (accounting for circular genome)
#             genome_length = len(self.protein_positions)
#             direct_distance = abs(pos1 - pos2)
#             circular_distance = genome_length - direct_distance
#             distance = min(direct_distance, circular_distance)
            
#             distances.append(distance)
#             pairs.append((prot1, prot2))
#             logFCmax_values.append(row['logFCmax'])
            
#         return {
#             'distances': distances,
#             'pairs': pairs,
#             'logFCmax': logFCmax_values
#         }
    
#     def generate_random_distances(self, n_samples: int = 1000) -> List[int]:
#         """Generate random protein pair distances for comparison"""
#         random_distances = []
#         proteins = list(self.protein_positions.keys())
        
#         for _ in range(n_samples):
#             # Randomly select two proteins
#             prot1, prot2 = np.random.choice(proteins, size=2, replace=False)
            
#             pos1 = self.protein_positions[prot1]
#             pos2 = self.protein_positions[prot2]
            
#             # Calculate distance same way as real pairs
#             genome_length = len(self.protein_positions)
#             direct_distance = abs(pos1 - pos2)
#             circular_distance = genome_length - direct_distance
#             distance = min(direct_distance, circular_distance)
            
#             random_distances.append(distance)
            
#         return random_distances
    
#     def analyze_proximity(self, n_random_samples: int = 1000) -> Dict:
#         """
#         Analyze if interacting proteins are closer in the genome than random pairs
        
#         Returns:
#         --------
#         dict : Analysis results including distances, statistics, and correlation
#         """
#         distance_data = self.calculate_genomic_distances()
#         random_distances = self.generate_random_distances(n_random_samples)
        
#         # Perform statistical test
#         statistic, p_value = stats.mannwhitneyu(
#             distance_data['distances'], 
#             random_distances,
#             alternative='less'  # Test if real distances are smaller than random
#         )
        
#         # Calculate correlation between distance and interaction strength
#         correlation, corr_p_value = stats.spearmanr(
#             distance_data['distances'],
#             distance_data['logFCmax']
#         )
        
#         results = {
#             'mean_distance': np.mean(distance_data['distances']),
#             'median_distance': np.median(distance_data['distances']),
#             'p_value': p_value,
#             'correlation_with_logFC': correlation,
#             'correlation_p_value': corr_p_value,
#             'n_interactions': len(distance_data['distances'])
#         }
        
#         return results
    
#     def plot_analysis(self, n_random_samples: int = 1000):
#         """Generate comprehensive visualization of the analysis"""
#         distance_data = self.calculate_genomic_distances()
#         random_distances = self.generate_random_distances(n_random_samples)
        
#         fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
        
#         # Plot 1: Distance distribution comparison
#         ax1.hist(random_distances, bins=30, alpha=0.5, label='Random pairs', density=True)
#         ax1.hist(distance_data['distances'], bins=30, alpha=0.5, label='Interacting pairs', density=True)
#         ax1.set_xlabel('Genomic distance')
#         ax1.set_ylabel('Density')
#         ax1.set_title('Distribution of genomic distances:\nInteracting vs Random protein pairs')
#         ax1.legend()
        
#         # Plot 2: Distance vs Interaction Strength
#         ax2.scatter(distance_data['distances'], distance_data['logFCmax'], alpha=0.5)
#         ax2.set_xlabel('Genomic distance')
#         ax2.set_ylabel('logFCmax (Interaction strength)')
#         ax2.set_title('Genomic Distance vs Interaction Strength')
        
#         plt.tight_layout()
#         plt.show()

# # Example usage:
# if __name__ == "__main__":
#     analyzer = ProteinProximityAnalyzer(
#         ppi_file="../Results/interactions_12Jan2025.xlsx",
#         genome_order_file="genome_order.csv"
#     )
    
#     # Run analysis
#     results = analyzer.analyze_proximity()
    
#     # Print results
#     print("\nAnalysis Results:")
#     print(f"Number of interactions analyzed: {results['n_interactions']}")
#     print(f"Mean genomic distance: {results['mean_distance']:.2f}")
#     print(f"Median genomic distance: {results['median_distance']:.2f}")
#     print(f"P-value (vs random): {results['p_value']:.2e}")
#     print(f"Correlation with interaction strength: {results['correlation_with_logFC']:.3f}")
#     print(f"Correlation p-value: {results['correlation_p_value']:.2e}")
    
#     # Generate visualizations
#     analyzer.plot_analysis()
########################################################################################
## ONLY ADJACENT
# import pandas as pd
# import numpy as np
# from scipy import stats
# from typing import List, Dict, Set, Tuple

# class ProteinAdjacencyAnalyzer:
#     def __init__(self, ppi_file: str, genome_order_file: str):
#         """
#         Initialize analyzer for adjacent protein interactions
#         """
#         self.ppis = self._load_interactions(ppi_file)
#         self.protein_order = self._load_genome_order(genome_order_file)
#         self.adjacent_pairs = self._find_adjacent_pairs()
        
#     def _load_interactions(self, ppi_file: str) -> List[Tuple[str, str]]:
#         """Load protein-protein interactions from Excel file"""
#         df = pd.read_excel(ppi_file)
#         self.full_ppi_data = df
#         return list(zip(df['Protein1'], df['Protein2']))
    
#     def _load_genome_order(self, genome_file: str) -> List[str]:
#         """
#         Load genomic order as a list of protein IDs in order of appearance
#         """
#         proteins_seen = set()
#         ordered_proteins = []
        
#         with open(genome_file, 'r') as f:
#             for line in f:
#                 if line.strip():
#                     protein = line.strip().split(':')[0]
#                     if protein not in proteins_seen:
#                         proteins_seen.add(protein)
#                         ordered_proteins.append(protein)
        
#         return ordered_proteins

#     def _find_adjacent_pairs(self) -> Set[Tuple[str, str]]:
#         """Find all adjacent protein pairs in genome"""
#         adjacent_pairs = set()
        
#         # Check each consecutive pair in the ordered list
#         for i in range(len(self.protein_order)-1):
#             p1, p2 = self.protein_order[i], self.protein_order[i+1]
#             adjacent_pairs.add((p1, p2))
#             adjacent_pairs.add((p2, p1))  # Add both orientations
                
#         return adjacent_pairs

#     def analyze_adjacency(self) -> Dict:
#         """
#         Analyze if interacting proteins are more likely to be adjacent
#         """
#         # Count interacting pairs that are adjacent
#         adjacent_interactions = sum(1 for pair in self.ppis 
#                                  if pair in self.adjacent_pairs)
        
#         total_interactions = len(self.ppis)
#         total_adjacent_pairs = len(self.adjacent_pairs) // 2  # Divide by 2 because we stored both orientations
        
#         # Calculate enrichment
#         fraction_adjacent = adjacent_interactions / total_interactions
#         expected_fraction = total_adjacent_pairs / (len(self.protein_order) * (len(self.protein_order) - 1) / 2)
#         enrichment = fraction_adjacent / expected_fraction if expected_fraction > 0 else 0
        
#         # Calculate statistics for adjacent vs non-adjacent interactions
#         adjacent_strengths = []
#         nonadjacent_strengths = []
        
#         adjacent_pairs_list = []
        
#         for idx, row in self.full_ppi_data.iterrows():
#             pair = (row['Protein1'], row['Protein2'])
#             if pair in self.adjacent_pairs:
#                 adjacent_strengths.append(row['logFCmax'])
#                 adjacent_pairs_list.append((row['Protein1'], row['Protein2'], row['logFCmax']))
#             else:
#                 nonadjacent_strengths.append(row['logFCmax'])
                
#         # Sort adjacent pairs by interaction strength
#         adjacent_pairs_list.sort(key=lambda x: x[2], reverse=True)
                
#         # Perform statistical test
#         if adjacent_strengths and nonadjacent_strengths:
#             statistic, pvalue = stats.mannwhitneyu(
#                 adjacent_strengths,
#                 nonadjacent_strengths,
#                 alternative='two-sided'
#             )
#         else:
#             statistic, pvalue = 0, 1.0
            
#         results = {
#             'adjacent_interactions': adjacent_interactions,
#             'total_interactions': total_interactions,
#             'fraction_adjacent': fraction_adjacent,
#             'enrichment': enrichment,
#             'adjacent_mean_strength': np.mean(adjacent_strengths) if adjacent_strengths else 0,
#             'nonadjacent_mean_strength': np.mean(nonadjacent_strengths) if nonadjacent_strengths else 0,
#             'strength_comparison_pvalue': pvalue,
#             'adjacent_pairs': adjacent_pairs_list
#         }
        
#         return results

# if __name__ == "__main__":
#     analyzer = ProteinAdjacencyAnalyzer(
#         ppi_file="../Results/interactions_12Jan2025.xlsx",
#         genome_order_file="genome_order.csv"
#     )
    
#     results = analyzer.analyze_adjacency()
    
#     print("\nAdjacency Analysis Results:")
#     print(f"Adjacent interactions: {results['adjacent_interactions']}")
#     print(f"Total interactions: {results['total_interactions']}")
#     print(f"Fraction adjacent: {results['fraction_adjacent']:.4f}")
#     print(f"Enrichment over random: {results['enrichment']:.2f}x")
#     print(f"\nInteraction Strength Analysis:")
#     print(f"Mean strength (adjacent): {results['adjacent_mean_strength']:.2f}")
#     print(f"Mean strength (non-adjacent): {results['nonadjacent_mean_strength']:.2f}")
#     print(f"Strength comparison p-value: {results['strength_comparison_pvalue']:.2e}")
    
#     print("\nTop adjacent interacting pairs by strength:")
#     for p1, p2, strength in results['adjacent_pairs'][:10]:
#         print(f"{p1} - {p2}: {strength:.2f}")
#################################################################
## 30 POSITIONS AWAY IN THE GENOME
import pandas as pd
import numpy as np
from collections import defaultdict
from scipy import stats

class ProximateProteinAnalyzer:
    def __init__(self, ppi_file: str, genome_order_file: str, max_distance: int = 10):
        self.max_distance = max_distance
        self.ppis = self._load_interactions(ppi_file)
        self.protein_positions = self._load_genome_order(genome_order_file)
        self.proximate_pairs = self._find_proximate_pairs()
        
    def _load_interactions(self, ppi_file: str):
        return pd.read_excel(ppi_file)
    
    def _load_genome_order(self, genome_file: str):
        with open(genome_file, 'r') as f:
            lines = f.readlines()
        
        positions = {}
        ordered_proteins = []
        for line in lines:
            if line.strip():
                protein = line.strip()
                if protein not in positions:
                    positions[protein] = len(ordered_proteins)
                    ordered_proteins.append(protein)
        
        return positions, ordered_proteins
    
    def _find_proximate_pairs(self):
        """Find all protein pairs within max_distance in genome"""
        positions, ordered_proteins = self.protein_positions
        genome_length = len(ordered_proteins)
        proximate_pairs = set()
        
        for i in range(genome_length):
            # Forward direction
            for j in range(1, self.max_distance + 1):
                pos2 = (i + j) % genome_length
                pair = (ordered_proteins[i], ordered_proteins[pos2])
                proximate_pairs.add(pair)
                proximate_pairs.add(pair[::-1])
                
            # Backward direction (for circular genome)
            for j in range(1, self.max_distance + 1):
                pos2 = (i - j) % genome_length
                pair = (ordered_proteins[i], ordered_proteins[pos2])
                proximate_pairs.add(pair)
                proximate_pairs.add(pair[::-1])
        
        return proximate_pairs
    
    def calculate_enrichment(self):
        """Calculate enrichment of proximate protein interactions"""
        positions, ordered_proteins = self.protein_positions
        genome_length = len(ordered_proteins)
        
        # Count observed proximate interactions
        proximate_interactions = set()
        for _, row in self.ppis.iterrows():
            prot1, prot2 = row['Protein1'], row['Protein2']
            if (prot1, prot2) in self.proximate_pairs:
                proximate_interactions.add((prot1, prot2))
        
        n_observed_proximate = len(proximate_interactions)
        total_interactions = len(self.ppis)
        
        # Calculate expected proportion
        total_possible_pairs = (genome_length * (genome_length - 1)) // 2
        total_possible_proximate = len(self.proximate_pairs) // 2  # Divide by 2 because pairs are bidirectional
        expected_proportion = total_possible_proximate / total_possible_pairs
        expected_proximate = total_interactions * expected_proportion
        
        # Calculate enrichment
        enrichment_fold = (n_observed_proximate / total_interactions) / expected_proportion if expected_proportion > 0 else 0
        
        # Statistical test using Fisher's exact test
        contingency_table = [
            [n_observed_proximate, total_interactions - n_observed_proximate],
            [total_possible_proximate, total_possible_pairs - total_possible_proximate]
        ]
        odds_ratio, p_value = stats.fisher_exact(contingency_table)
        
        results = {
            'n_observed_proximate': n_observed_proximate,
            'total_interactions': total_interactions,
            'expected_proximate': expected_proximate,
            'enrichment_fold': enrichment_fold,
            'p_value': p_value,
            'odds_ratio': odds_ratio,
            'total_possible_proximate': total_possible_proximate,
            'total_possible_pairs': total_possible_pairs
        }

        print("\nEnrichment Analysis Results:")
        print(f"Observed proximate interactions: {results['n_observed_proximate']}")
        print(f"Expected proximate interactions: {results['expected_proximate']:.1f}")
        print(f"Enrichment fold: {results['enrichment_fold']:.2f}")
        print(f"Fisher's exact test p-value: {results['p_value']:.2e}")
        print(f"Odds ratio: {results['odds_ratio']:.2f}")
        
        return results
    
    def analyze_proximity_vs_random(self, n_permutations=1000):
        """Compare proximate interactions to random expectation using permutation test"""
        positions, ordered_proteins = self.protein_positions
        
        # Observed number of proximate interactions
        observed_proximate = sum(1 for _, row in self.ppis.iterrows() 
                               if (row['Protein1'], row['Protein2']) in self.proximate_pairs)
        
        # Permutation test
        random_counts = []
        proteins = list(positions.keys())
        
        for _ in range(n_permutations):
            # Shuffle protein positions
            shuffled_positions = dict(zip(proteins, np.random.permutation(list(positions.values()))))
            
            # Count proximate interactions with shuffled positions
            proximate_count = 0
            for _, row in self.ppis.iterrows():
                prot1, prot2 = row['Protein1'], row['Protein2']
                if prot1 in shuffled_positions and prot2 in shuffled_positions:
                    pos1 = shuffled_positions[prot1]
                    pos2 = shuffled_positions[prot2]
                    
                    # Calculate distance considering circular genome
                    direct_distance = abs(pos1 - pos2)
                    circular_distance = len(proteins) - direct_distance
                    distance = min(direct_distance, circular_distance)
                    
                    if distance <= self.max_distance:
                        proximate_count += 1
            
            random_counts.append(proximate_count)
        
        # Calculate p-value
        random_counts = np.array(random_counts)
        p_value = np.mean(random_counts >= observed_proximate)
        
        z_score = (observed_proximate - np.mean(random_counts)) / np.std(random_counts)
        
        permutation_results = {
            'observed_proximate': observed_proximate,
            'random_mean': np.mean(random_counts),
            'random_std': np.std(random_counts),
            'z_score': z_score,
            'p_value': p_value,
            'n_permutations': n_permutations
        }
        
        return permutation_results

    def print_enrichment_results(self, enrichment_results, permutation_results):
        """Print analysis results"""
        print("\nEnrichment Analysis Results:")
        print(f"Observed proximate interactions: {enrichment_results['n_observed_proximate']}")
        print(f"Expected proximate interactions: {enrichment_results['expected_proximate']:.1f}")
        print(f"Enrichment fold: {enrichment_results['enrichment_fold']:.2f}")
        print(f"Fisher's exact test p-value: {enrichment_results['p_value']:.2e}")
        print(f"Odds ratio: {enrichment_results['odds_ratio']:.2f}")
        
        print("\nPermutation Test Results:")
        print(f"Observed proximate interactions: {permutation_results['observed_proximate']}")
        print(f"Random expectation: {permutation_results['random_mean']:.1f} ± {permutation_results['random_std']:.1f}")
        print(f"Z-score: {permutation_results['z_score']:.2f}")
        print(f"Permutation test p-value: {permutation_results['p_value']:.2e}")

if __name__ == "__main__":
    analyzer = ProximateProteinAnalyzer(
        ppi_file="../Results/protProt_Libs2_logFC1.5FDR0.05_woutSticky_FragPairs2.xlsx",
        genome_order_file="genome_order.csv",
        max_distance=3
    )
    
    enrichment_results = analyzer.calculate_enrichment()
    print("Calculated enrichment....")
    permutation_results = analyzer.analyze_proximity_vs_random()
    print("Proximity vs random done....")
    analyzer.print_enrichment_results(enrichment_results, permutation_results)