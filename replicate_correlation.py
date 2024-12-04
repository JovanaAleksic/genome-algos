import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

def load_and_process_data(file1, file2):
    """Load and process the two replicate files"""
    df1 = pd.read_csv(file1, sep='\t')
    df2 = pd.read_csv(file2, sep='\t')
    
    # Make sure we're using the genes as index
    df1 = df1.set_index('genes')
    df2 = df2.set_index('genes')
    
    return df1, df2

def compare_replicates(df1, df2, fdr_threshold):
    """Compare replicates at given FDR threshold"""
    # Filter based on FDR threshold for both datasets
    df1_filtered = df1[df1['FDR'] < fdr_threshold]
    df2_filtered = df2[df2['FDR'] < fdr_threshold]
    
    # Find common protein pairs
    common_pairs = df1_filtered.index.intersection(df2_filtered.index)
    
    # Get logFC values for common pairs
    logFC1 = df1_filtered.loc[common_pairs, 'logFC']
    logFC2 = df2_filtered.loc[common_pairs, 'logFC']
    
    # Calculate R-squared
    slope, intercept, r_value, p_value, std_err = stats.linregress(logFC1, logFC2)
    r_squared = r_value**2
    
    return logFC1, logFC2, r_squared, common_pairs

def create_correlation_plot(fdr_thresholds):
    """Create correlation plots for different FDR thresholds"""
    df1, df2 = load_and_process_data('../Counts/Diff/20230707_FLPRS_LS_1.2mMATdiff.txt', 
                                    '../Counts/Diff/20230707_FLPRS_LS_2.2mMATdiff.txt')
    
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    fig.suptitle('Correlation of logFC Values Between Replicates at Different FDR Thresholds')
    
    for i, fdr in enumerate(fdr_thresholds):
        logFC1, logFC2, r_squared, common_pairs = compare_replicates(df1, df2, fdr)
        
        # Create scatter plot
        axes[i].scatter(logFC1, logFC2, alpha=0.5)
        
        # Add diagonal line
        min_val = min(min(logFC1), min(logFC2))
        max_val = max(max(logFC1), max(logFC2))
        axes[i].plot([min_val, max_val], [min_val, max_val], 'r--', alpha=0.5)
        
        # Add labels and title
        axes[i].set_xlabel('Replicate 1 logFC')
        axes[i].set_ylabel('Replicate 2 logFC')
        axes[i].set_title(f'FDR < {fdr}\nR² = {r_squared:.3f}')
        axes[i].grid(True, alpha=0.3)
        
        # Make plot square
        axes[i].set_aspect('equal', adjustable='box')
    
    plt.tight_layout()
    return plt

# Create plots for specified FDR thresholds
fdr_thresholds = [0.1, 0.2, 0.4]
plot = create_correlation_plot(fdr_thresholds)
plot.show()

# Print detailed statistics
df1, df2 = load_and_process_data('../Counts/Diff/20230707_FLPRS_LS_1.2mMATdiff.txt', 
                                '../Counts/Diff/20230707_FLPRS_LS_2.2mMATdiff.txt')