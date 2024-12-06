# Matthew Beltran
# 11/25/2024
# Using chatgpt code copilot

import MAGE
import numpy as np

def test():
    # --- Generate synthetic data for testing ---
    num_genes = 100
    num_replicates = 10
    corr = 0.9

    # Random gene expression data for X and Y (simulating biological data)
    #data_x = np.random.rand(num_genes, num_replicates) * 100
    #data_y = np.random.rand(num_genes, num_replicates) * 100
    data_x, data_y = generate_data(num_genes, num_replicates, corr)
    # --- Run MAGE function ---
    print("Running MAGE...")
    density_matrix = MAGE.mage(data_x, data_y)

def generate_data(num_genes, replicates, correlation=0.5, avg_std_x=5, avg_std_y=5):
    """
    Generates synthetic data_x and data_y with a specified positive correlation 
    and adjustable average standard deviations (STD) for genes.

    Parameters:
    - num_genes: int
        Number of genes.
    - replicates: int
        Number of replicates per gene.
    - correlation: float
        Desired correlation coefficient between data_x and data_y.
    - avg_std_x: float
        Desired average standard deviation for data_x.
    - avg_std_y: float
        Desired average standard deviation for data_y.

    Returns:
    - data_x: np.ndarray
        Simulated gene expression data for data_x.
    - data_y: np.ndarray
        Simulated gene expression data for data_y.
    """
    # Mean and covariance for multivariate normal distribution
    mean = [0, 0]  # Mean of X and Y
    covariance = [[1, correlation], [correlation, 1]]  # Covariance matrix

    # Generate correlated data
    raw_data = np.random.multivariate_normal(mean, covariance, num_genes * replicates)

    # Reshape into (genes x replicates)
    data_x = raw_data[:, 0].reshape(num_genes, replicates)
    data_y = raw_data[:, 1].reshape(num_genes, replicates)

    # Scale to a positive range for expression values
    data_x = (data_x - data_x.min()) / (data_x.max() - data_x.min()) * 100
    data_y = (data_y - data_y.min()) / (data_y.max() - data_y.min()) * 100

    # Adjust the standard deviation of each gene in data_x
    current_avg_std_x = np.mean(np.std(data_x, axis=1))
    scale_factor_x = avg_std_x / current_avg_std_x
    data_x = np.mean(data_x, axis=1, keepdims=True) + scale_factor_x * (data_x - np.mean(data_x, axis=1, keepdims=True))

    # Adjust the standard deviation of each gene in data_y
    current_avg_std_y = np.mean(np.std(data_y, axis=1))
    scale_factor_y = avg_std_y / current_avg_std_y
    data_y = np.mean(data_y, axis=1, keepdims=True) + scale_factor_y * (data_y - np.mean(data_y, axis=1, keepdims=True))

    return data_x, data_y

if __name__ == "__main__":
    test()