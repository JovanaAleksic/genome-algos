import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


def load_data_files(file_paths):
    """Load all replicate files and return as a list of dataframes"""
    dataframes = []
    for file in file_paths:
        df = pd.read_csv(file, sep='\t')
        dataframes.append(df)
    return dataframes

def find_common_top_interactions(dataframes, n=50):
    """Find common top interactions across all replicates based on logFC"""
    # Get logFC and sort for each dataframe
    top_pairs = []
    for df in dataframes:
        sorted_df = df.nlargest(300, 'logFC')  # Get more than 10 to ensure overlap
        top_pairs.append(set(sorted_df['genes']))
    
    # Find common pairs across all replicates
    common_pairs = set.intersection(*map(set, top_pairs))
    
    # Get mean absolute logFC for common pairs to rank them
    mean_logFC = {}
    for pair in common_pairs:
        mean_abs_fc = np.mean([df.loc[df['genes'] == pair, 'logFC'].iloc[0]
                              for df in dataframes])
        mean_logFC[pair] = mean_abs_fc
    
    # Get top N common pairs
    top_common = sorted(mean_logFC.items(), key=lambda x: x[1], reverse=True)[:n]
    indexes = [0,1,6,9,10,24,28,49]
    top_common = [top_common[i] for i in indexes]
    return [pair[0] for pair in top_common]


def create_volcano_plots(file_paths, fc_cutoff=0.8, p_cutoff=0.1, y_max=20):
    """Create three volcano plots with common top interactions highlighted"""
    
    # Load data
    dataframes = load_data_files(file_paths)
    
    # Find common top interactions
    top_common = find_common_top_interactions(dataframes)
    
    # Define default label positions for each plot
    # Format: {replicate_number: {gene_pair: (x_offset, y_offset)}}
    label_positions = {
        0: {  # First plot
            top_common[0]: (-60, -15),
            top_common[2]: (0, -15),
            top_common[1]: (-85,-8),
            top_common[3]: (-15, 7),
            top_common[4]: (-10, 10),
            top_common[5]: (10, -3),
            top_common[6]: (-10, -12),
            top_common[7]: (-15, -15)
        },
        1: {  # Second plot
            top_common[0]: (-50, -15),
            top_common[2]: (0, -15),
            top_common[1]: (-90, -8),
            top_common[3]: (-15, 7),
            top_common[4]: (-10, -12),
            top_common[5]: (-10, -12),
            top_common[6]: (-10, 7),
            top_common[7]: (-15, -15)
        },
        2: {  # Third plot
            top_common[0]: (15, -15),
            top_common[2]: (0, -15),
            top_common[1]: (-95, -10),
            top_common[3]: (-15, 7),
            top_common[4]: (-10, 10),
            top_common[5]: (10, -3),
            top_common[6]: (10,0),
            top_common[7]: (-15, -15)
        }
    }
    
    # Create figure with three subplots
    fig, axes = plt.subplots(1, 3, figsize=(20, 7))
    fig.suptitle('Volcano Plots of Protein Pair Interactions Across Replicates', fontsize=14)
    
    # Create volcano plot for each replicate
    for i, (df, ax) in enumerate(zip(dataframes, axes)):
        # Calculate -log10(FDR) and cap at y_max
        df['neg_log_fdr'] = -np.log10(df['FDR'])
        df['neg_log_fdr'] = df['neg_log_fdr'].clip(upper=y_max)
        
        # Define significance thresholds
        significant = (df['logFC'] >= fc_cutoff) & (df['FDR'] < p_cutoff)

        df = df[df['logFC'] >= 0]
        
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
                ax.scatter(point['logFC'], min(point['neg_log_fdr'], y_max),
                         c='gold', s=100, edgecolor='black', zorder=5)
                
                # Get label position for this pair in this plot
                x_offset, y_offset = label_positions[i][pair]
                
                ax.annotate(pair.split(':')[0] + ':' + pair.split(':')[1],
                          (point['logFC'], min(point['neg_log_fdr'], y_max)),
                          xytext=(x_offset, y_offset), 
                          textcoords='offset points',
                          fontsize=8, fontweight='bold', color='black')
        
        # Add lines for cutoffs
        ax.axhline(y=-np.log10(p_cutoff), color='blue', linestyle='--', alpha=0.5)
        ax.axvline(x=fc_cutoff, color='blue', linestyle='--', alpha=0.5)
        
        # Customize the plot
        ax.set_xlabel('Log2 Fold Change')
        ax.set_ylabel('-log10(FDR)' if i == 0 else '')
        ax.set_title(f'Replicate {i+1}')
        ax.set_ylim(0, y_max)
        ax.grid(True, alpha=0.3)
        
        # Add legend to first plot only
        if i == 0:
            handles, labels = ax.get_legend_handles_labels()
            handles.append(plt.scatter([], [], c='gold', s=100, edgecolor='black'))
            handles.append(plt.scatter([], [], marker='^', c='black', s=50))
            labels.append('Highlighted Interactions')
            ax.legend(handles, labels, loc='lower right')
    
    plt.tight_layout()
    # Save as SVG
    plt.savefig('volcano_oneSided.svg', format='svg', dpi=300, bbox_inches='tight')

    return plt, top_common

# File paths
file_paths = [
    '../Counts/Diff/20230707_FLPRS_LS_1.2mMATdiff.txt',
    '../Counts/Diff/20230707_FLPRS_LS_2.2mMATdiff.txt',
    '../Counts/Diff/20230707_FLPRS_LS_3.2mMATdiff.txt'
]

# Create plots and get top common interactions
plot, top_common = create_volcano_plots(file_paths)

# Print top common interactions
print("\nTop 10 common interactions across all replicates:")
print("-" * 50)
for i, pair in enumerate(top_common, 1):
    print(f"{i}. {pair}")


plot.show()


