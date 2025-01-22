import pandas as pd
import numpy as np
from scipy import stats
from typing import List, Tuple, Dict
import matplotlib.pyplot as plt

class ProteinProximityAnalyzer:
    def __init__(self, ppi_file: str, genome_order_file: str):
        """
        Initialize the analyzer with protein interaction and genome order data
        
        Parameters:
        -----------
        ppi_file : str
            Path to Excel file containing protein-protein interactions
            Expected columns: Protein1, Protein2, plus additional metadata
        genome_order_file : str
            Path to file containing protein order in genome
            Expected format: one column with protein IDs in genomic order
        """
        self.ppis = self._load_interactions(ppi_file)
        self.protein_positions = self._load_genome_order(genome_order_file)
        
    def _load_interactions(self, ppi_file: str) -> List[Tuple[str, str]]:
        """Load protein-protein interactions from Excel file"""
        # Read Excel file
        df = pd.read_excel(ppi_file)
        
        # Extract just the protein pairs
        protein_pairs = list(zip(df['Protein1'], df['Protein2']))
        
        # Store the full dataframe for later analysis if needed
        self.full_ppi_data = df
        
        return protein_pairs
    
    def _load_genome_order(self, genome_file: str) -> Dict[str, int]:
        """Load genomic positions and create position lookup dictionary"""
        df = pd.read_csv(genome_file)
        # Assuming single column with protein IDs
        proteins = df.iloc[:, 0].tolist()
        return {protein: pos for pos, protein in enumerate(proteins)}
    
    def calculate_genomic_distances(self) -> Dict[str, List]:
        """Calculate genomic distances between all interacting proteins"""
        distances = []
        pairs = []
        logFCmax_values = []
        
        for idx, row in self.full_ppi_data.iterrows():
            prot1, prot2 = row['Protein1'], row['Protein2']
            
            # Skip if either protein is not in the genome order list
            if prot1 not in self.protein_positions or prot2 not in self.protein_positions:
                continue
                
            pos1 = self.protein_positions[prot1]
            pos2 = self.protein_positions[prot2]
            
            # Calculate shortest distance (accounting for circular genome)
            genome_length = len(self.protein_positions)
            direct_distance = abs(pos1 - pos2)
            circular_distance = genome_length - direct_distance
            distance = min(direct_distance, circular_distance)
            
            distances.append(distance)
            pairs.append((prot1, prot2))
            logFCmax_values.append(row['logFCmax'])
            
        return {
            'distances': distances,
            'pairs': pairs,
            'logFCmax': logFCmax_values
        }
    
    def generate_random_distances(self, n_samples: int = 1000) -> List[int]:
        """Generate random protein pair distances for comparison"""
        random_distances = []
        proteins = list(self.protein_positions.keys())
        
        for _ in range(n_samples):
            # Randomly select two proteins
            prot1, prot2 = np.random.choice(proteins, size=2, replace=False)
            
            pos1 = self.protein_positions[prot1]
            pos2 = self.protein_positions[prot2]
            
            # Calculate distance same way as real pairs
            genome_length = len(self.protein_positions)
            direct_distance = abs(pos1 - pos2)
            circular_distance = genome_length - direct_distance
            distance = min(direct_distance, circular_distance)
            
            random_distances.append(distance)
            
        return random_distances
    
    def analyze_proximity(self, n_random_samples: int = 1000) -> Dict:
        """
        Analyze if interacting proteins are closer in the genome than random pairs
        
        Returns:
        --------
        dict : Analysis results including distances, statistics, and correlation
        """
        distance_data = self.calculate_genomic_distances()
        random_distances = self.generate_random_distances(n_random_samples)
        
        # Perform statistical test
        statistic, p_value = stats.mannwhitneyu(
            distance_data['distances'], 
            random_distances,
            alternative='less'  # Test if real distances are smaller than random
        )
        
        # Calculate correlation between distance and interaction strength
        correlation, corr_p_value = stats.spearmanr(
            distance_data['distances'],
            distance_data['logFCmax']
        )
        
        results = {
            'mean_distance': np.mean(distance_data['distances']),
            'median_distance': np.median(distance_data['distances']),
            'p_value': p_value,
            'correlation_with_logFC': correlation,
            'correlation_p_value': corr_p_value,
            'n_interactions': len(distance_data['distances'])
        }
        
        return results
    
    def plot_analysis(self, n_random_samples: int = 1000):
        """Generate comprehensive visualization of the analysis"""
        distance_data = self.calculate_genomic_distances()
        random_distances = self.generate_random_distances(n_random_samples)
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
        
        # Plot 1: Distance distribution comparison
        ax1.hist(random_distances, bins=30, alpha=0.5, label='Random pairs', density=True)
        ax1.hist(distance_data['distances'], bins=30, alpha=0.5, label='Interacting pairs', density=True)
        ax1.set_xlabel('Genomic distance')
        ax1.set_ylabel('Density')
        ax1.set_title('Distribution of genomic distances:\nInteracting vs Random protein pairs')
        ax1.legend()
        
        # Plot 2: Distance vs Interaction Strength
        ax2.scatter(distance_data['distances'], distance_data['logFCmax'], alpha=0.5)
        ax2.set_xlabel('Genomic distance')
        ax2.set_ylabel('logFCmax (Interaction strength)')
        ax2.set_title('Genomic Distance vs Interaction Strength')
        
        plt.tight_layout()
        plt.show()

# Example usage:
if __name__ == "__main__":
    analyzer = ProteinProximityAnalyzer(
        ppi_file="../Results/interactions_12Jan2025.xlsx",
        genome_order_file="genome_order.csv"
    )
    
    # Run analysis
    results = analyzer.analyze_proximity()
    
    # Print results
    print("\nAnalysis Results:")
    print(f"Number of interactions analyzed: {results['n_interactions']}")
    print(f"Mean genomic distance: {results['mean_distance']:.2f}")
    print(f"Median genomic distance: {results['median_distance']:.2f}")
    print(f"P-value (vs random): {results['p_value']:.2e}")
    print(f"Correlation with interaction strength: {results['correlation_with_logFC']:.3f}")
    print(f"Correlation p-value: {results['correlation_p_value']:.2e}")
    
    # Generate visualizations
    analyzer.plot_analysis()