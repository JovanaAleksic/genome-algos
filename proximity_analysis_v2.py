import pandas as pd
import numpy as np
from collections import defaultdict
from scipy import stats

class InteractionNeighborAnalyzer:
    def __init__(self, ppi_file: str, genome_order_file: str, max_distance: int = 3):
        self.max_distance = max_distance
        self.ppis = self._load_interactions(ppi_file)
        self.protein_positions, self.ordered_proteins = self._load_genome_order(genome_order_file)
        self.interaction_partners = self._get_interaction_partners()
        
    def _load_interactions(self, ppi_file: str):
        return pd.read_excel(ppi_file)
    
    def _load_genome_order(self, genome_file: str):
        """Load protein positions from genome order file"""
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
    
    def _get_interaction_partners(self):
        """Create dictionary of interaction partners for each protein"""
        partners = defaultdict(set)
        for _, row in self.ppis.iterrows():
            prot1, prot2 = row['Fragment1'], row['Fragment2']
            partners[prot1].add(prot2)
            partners[prot2].add(prot1)
        return partners
    
    def _get_genomic_neighbors(self, protein):
        """Get proteins within max_distance in genome"""
        if protein not in self.protein_positions:
            return set()
            
        pos = self.protein_positions[protein]
        genome_length = len(self.ordered_proteins)
        neighbors = set()
        
        # Forward direction
        for i in range(1, self.max_distance + 1):
            neighbor_pos = (pos + i) % genome_length
            neighbors.add(self.ordered_proteins[neighbor_pos])
            
        # Backward direction
        for i in range(1, self.max_distance + 1):
            neighbor_pos = (pos - i) % genome_length
            neighbors.add(self.ordered_proteins[neighbor_pos])
            
        return neighbors
    
    def analyze_neighbor_interactions(self):
        """Analyze enrichment of interactions with partners' genomic neighbors"""
        total_observed = 0
        total_possible = 0
        
        # For each protein and its interaction partners
        for protein in self.interaction_partners:
            # Get protein's direct interaction partners
            direct_partners = self.interaction_partners[protein]
            
            # For each partner, get its genomic neighbors and check if protein interacts with them
            for partner in direct_partners:
                # Get genomic neighbors of this partner
                partner_neighbors = self._get_genomic_neighbors(partner)
                
                # Remove the original protein and the partner itself
                partner_neighbors -= {protein, partner}
                
                # Get the protein's interaction partners that are also neighbors of its partner
                neighbor_interactions = len(direct_partners.intersection(partner_neighbors))
                
                total_observed += neighbor_interactions
                total_possible += len(partner_neighbors)
        
        # Calculate enrichment statistics
        if total_possible == 0:
            return {
                'enrichment': 0,
                'total_observed': 0,
                'total_possible': 0,
                'background_rate': 0,
                'p_value': 1.0
            }
            
        # Calculate background interaction rate
        total_proteins = len(self.protein_positions)
        total_interactions = sum(len(partners) for partners in self.interaction_partners.values()) / 2
        background_rate = total_interactions / (total_proteins * (total_proteins - 1) / 2)
        
        # Expected number of interactions given background rate
        expected = background_rate * total_possible
        
        # Calculate enrichment
        enrichment = (total_observed / total_possible) / background_rate if background_rate > 0 else 0
        
        # Statistical test using Fisher's exact test
        contingency_table = [
            [total_observed, total_possible - total_observed],
            [total_interactions, (total_proteins * (total_proteins - 1) / 2) - total_interactions]
        ]
        _, p_value = stats.fisher_exact(contingency_table)
        
        return {
            'enrichment': enrichment,
            'total_observed': total_observed,
            'total_possible': total_possible,
            'expected': expected,
            'background_rate': background_rate,
            'p_value': p_value
        }
    
    def print_results(self, results):
        """Print analysis results"""
        print("\nInteraction Partner Neighbor Analysis Results:")
        print(f"Total observed neighbor interactions: {results['total_observed']}")
        print(f"Total possible neighbor interactions: {results['total_possible']}")
        print(f"Expected neighbor interactions: {results['expected']:.1f}")
        print(f"Background interaction rate: {results['background_rate']:.3f}")
        print(f"Enrichment fold: {results['enrichment']:.2f}")
        print(f"Fisher's exact test p-value: {results['p_value']:.2e}")

# Example usage:
if __name__ == "__main__":
    analyzer = InteractionNeighborAnalyzer(
        ppi_file="../Results/fragFrag_logFC1.5FDR0.05.xlsx",
        genome_order_file="sorted_genomeFrags.txt",
        max_distance=3
    )
    
    results = analyzer.analyze_neighbor_interactions()
    analyzer.print_results(results)