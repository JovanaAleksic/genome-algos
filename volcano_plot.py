import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def create_volcano_plot(data, fc_cutoff=0.8, p_cutoff=0.1):
    """
    Create a volcano plot from differential analysis data
    
    """
    
    # Calculate -log10(FDR)
    data['neg_log_fdr'] = -np.log10(data['FDR'])
    
    # Create the scatter plot
    plt.figure(figsize=(12, 8))
    
    # Define significance thresholds
    significant = (data['logFC'] >= fc_cutoff) & (data['FDR'] < p_cutoff)
    
    # Plot non-significant points
    plt.scatter(data.loc[~significant, 'logFC'], 
               data.loc[~significant, 'neg_log_fdr'],
               c='grey', alpha=0.5, s=20, label='No Interaction')
    
    # Plot significant points
    plt.scatter(data.loc[significant, 'logFC'],
               data.loc[significant, 'neg_log_fdr'],
               c='red', alpha=0.7, s=20, label='Interaction')
    
    # Add horizontal line for p-value cutoff
    plt.axhline(y=-np.log10(p_cutoff), color='blue', linestyle='--', alpha=0.5)
    
    # Add vertical lines for fold change cutoff
    plt.axvline(x=fc_cutoff, color='blue', linestyle='--', alpha=0.5)
    
    # Customize the plot
    plt.xlabel('Log2 Fold Change')
    plt.ylabel('-log10(FDR)')
    plt.title('Volcano Plot of Protein Pair Interactions')
    
    # Add grid

    plt.grid(True, alpha=0.3)
    plt.legend(loc='upper right')

    plt.ylim(-1, 20)
    
    # top_significant = data[significant].nlargest(10, 'neg_log_fdr')
    # for idx, row in top_significant.iterrows():
    #     plt.annotate(row.name, 
    #                 (row['logFC'], row['neg_log_fdr']),
    #                 xytext=(5, 5), textcoords='offset points',
    #                 fontsize=8, alpha=0.7)
    
    plt.tight_layout()
    return plt

# Read the data
data = pd.read_csv('Counts/Diff/20230707_FLPRS_LS_1.2mMATdiff.txt', sep='\t')

# Create and show the plot
volcano_plot = create_volcano_plot(data)
volcano_plot.show()