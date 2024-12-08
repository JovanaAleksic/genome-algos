import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# def create_volcano_plot(data, fc_cutoff=0.8, p_cutoff=0.1):
#     """
#     Create a volcano plot from differential analysis data
    
#     """
    
#     # Calculate -log10(FDR)
#     data['neg_log_fdr'] = -np.log10(data['FDR'])
    
#     # Create the scatter plot
#     plt.figure(figsize=(12, 8))
    
#     # Define significance thresholds
#     significant = (data['logFC'] >= fc_cutoff) & (data['FDR'] < p_cutoff)
    
#     # Plot non-significant points
#     plt.scatter(data.loc[~significant, 'logFC'], 
#                data.loc[~significant, 'neg_log_fdr'],
#                c='grey', alpha=0.5, s=20, label='No Interaction')
    
#     # Plot significant points
#     plt.scatter(data.loc[significant, 'logFC'],
#                data.loc[significant, 'neg_log_fdr'],
#                c='red', alpha=0.7, s=20, label='Interaction')
    
#     # Add horizontal line for p-value cutoff
#     plt.axhline(y=-np.log10(p_cutoff), color='blue', linestyle='--', alpha=0.5)
    
#     # Add vertical lines for fold change cutoff
#     plt.axvline(x=fc_cutoff, color='blue', linestyle='--', alpha=0.5)
    
#     # Customize the plot
#     plt.xlabel('Log2 Fold Change')
#     plt.ylabel('-log10(FDR)')
#     plt.title('Volcano Plot of Protein Pair Interactions')
    
#     # Add grid

#     plt.grid(True, alpha=0.3)
#     plt.legend(loc='upper right')

#     plt.ylim(-1, 20)
    
#     # top_significant = data[significant].nlargest(10, 'neg_log_fdr')
#     # for idx, row in top_significant.iterrows():
#     #     plt.annotate(row.name, 
#     #                 (row['logFC'], row['neg_log_fdr']),
#     #                 xytext=(5, 5), textcoords='offset points',
#     #                 fontsize=8, alpha=0.7)
    
#     plt.tight_layout()
#     return plt

# # Read the data
# data = pd.read_csv('Counts/Diff/20230707_FLPRS_LS_1.2mMATdiff.txt', sep='\t')

# # Create and show the plot
# volcano_plot = create_volcano_plot(data)
# volcano_plot.show()



def load_data_files(file_paths):
    """Load all replicate files and return as a list of dataframes"""
    dataframes = []
    for file in file_paths:
        df = pd.read_csv(file, sep='\t')
        dataframes.append(df)
    return dataframes

def find_common_top_interactions(dataframes, n=10):
    """Find common top interactions across all replicates based on absolute logFC"""
    # Get absolute logFC and sort for each dataframe
    top_pairs = []
    for df in dataframes:
        sorted_df = df.nlargest(20, 'logFC')  # Get more than 10 to ensure overlap
        top_pairs.append(set(sorted_df['genes']))
    
    # Find common pairs across all replicates
    common_pairs = set.intersection(*map(set, top_pairs))
    
    # Get mean absolute logFC for common pairs to rank them
    mean_logFC = {}
    for pair in common_pairs:
        mean_abs_fc = np.mean([abs(df.loc[df['genes'] == pair, 'logFC'].iloc[0]) 
                              for df in dataframes])
        mean_logFC[pair] = mean_abs_fc
    
    # Get top N common pairs
    top_common = sorted(mean_logFC.items(), key=lambda x: x[1], reverse=True)[:n]
    return [pair[0] for pair in top_common]

def create_volcano_plots(file_paths, fc_cutoff=0.8, p_cutoff=0.1):
    """Create three volcano plots with common top interactions highlighted"""
    
    # Load data
    dataframes = load_data_files(file_paths)
    
    # Find common top interactions
    top_common = find_common_top_interactions(dataframes)
    
    # Create figure with three subplots
    fig, axes = plt.subplots(1, 3, figsize=(20, 7))
    fig.suptitle('Volcano Plots of Protein Pair Interactions Across Replicates', fontsize=14)
    
    # Create volcano plot for each replicate
    for i, (df, ax) in enumerate(zip(dataframes, axes)):
        # Calculate -log10(FDR)
        df['neg_log_fdr'] = -np.log10(df['FDR'])
        
        # Define significance thresholds
        significant = (df['logFC'] >= fc_cutoff) & (df['FDR'] < p_cutoff)
        
        # Plot non-significant points
        ax.scatter(df.loc[~significant, 'logFC'], 
                  df.loc[~significant, 'neg_log_fdr'],
                  c='grey', alpha=0.5, s=20, label='No Interaction')
        
        # Plot significant points
        ax.scatter(df.loc[significant, 'logFC'],
                  df.loc[significant, 'neg_log_fdr'],
                  c='red', alpha=0.7, s=20, label='Interaction')
        
        # Highlight common top interactions
        for pair in top_common:
            if pair in df['genes'].values:
                point = df[df['genes'] == pair].iloc[0]
                ax.scatter(point['logFC'], point['neg_log_fdr'],
                         c='gold', s=100, edgecolor='black', zorder=5)
                ax.annotate(pair.split(':')[0] + ':' + pair.split(':')[1][:3],
                          (point['logFC'], point['neg_log_fdr']),
                          xytext=(5, 5), textcoords='offset points',
                          fontsize=8, fontweight='bold', color='black')
        
        # Add lines for cutoffs
        ax.axhline(y=-np.log10(p_cutoff), color='blue', linestyle='--', alpha=0.5)
        ax.axvline(x=fc_cutoff, color='blue', linestyle='--', alpha=0.5)
        
        # Customize the plot
        ax.set_xlabel('Log2 Fold Change')
        ax.set_ylabel('-log10(FDR)' if i == 0 else '')
        ax.set_title(f'Replicate {i+1}')
        ax.set_ylim(0, 100)
        ax.grid(True, alpha=0.3)
        
        # Add legend to first plot only
        if i == 0:
            handles, labels = ax.get_legend_handles_labels()
            handles.append(plt.scatter([], [], c='gold', s=100, edgecolor='black'))
            labels.append('Common Top Interactions')
            ax.legend(handles, labels, loc='upper right')
    
    plt.tight_layout()
    return plt, top_common

# File paths
file_paths = [
    '../../Counts/Diff/20230707_FLPRS_LS_1.2mMATdiff.txt',
    '../../Counts/Diff/20230707_FLPRS_LS_2.2mMATdiff.txt',
    '../../Counts/Diff/20230707_FLPRS_LS_3.2mMATdiff.txt'
]

# Create plots and get top common interactions
plot, top_common = create_volcano_plots(file_paths)

# Print top common interactions
print("\nTop 10 common interactions across all replicates:")
print("-" * 50)
for i, pair in enumerate(top_common, 1):
    print(f"{i}. {pair}")

plot.show()