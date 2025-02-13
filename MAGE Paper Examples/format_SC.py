"""
    Reads an NCBI single-cell gene expression profile from a text file.
    
    Parameters:
        file_path (str): Path to the NCBI dataset file.
    
    Returns:
        sample_ids (list): List of sample IDs.
        gene_names (list): List of gene names.
        profile (ndarray): NumPy array containing gene expression data.
"""

import numpy as np
import pickle

file_path = r"C:\Users\Matthew\Documents\GitHub\MAGE\MAGE Paper Examples\workspaces\GSE98816_Brain_samples_normalized_counts_matrix.txt"
save_name = 'sc_brain.csv'

# Read the first line to extract sample IDs
with open(file_path, 'r') as file:
    first_line = file.readline().strip().split()
    sample_ids = [s.strip('"') for s in first_line[1:]]  # Remove quotes
    
# Count number of genes
num_genes = sum(1 for _ in open(file_path)) - 1  # Subtract header line
    
# Initialize storage
gene_names = []
profile = np.zeros((num_genes, len(sample_ids)))
    
# Read file line-by-line
with open(file_path, 'r') as file:
    next(file)  # Skip header
    for i, line in enumerate(file):
        parts = line.strip().split()
        gene_names.append(parts[0].strip('"'))  # Remove quotes
        profile[i, :] = np.array(parts[1:], dtype=float)  # Convert expression data to float
        print(f"Processing gene {i+1}/{num_genes}")
    
#return sample_ids, gene_names, profile
#sample_ids.insert(0,'gene')
#np.savetxt(save_name, np.c_[np.array(gene_names),profile], fmt='%s', delimiter=',',header=sample_ids)

with open("sc_brain.pkl", "wb") as f:
    pickle.dump((profile, gene_names, sample_ids), f)