# Matthew Beltran
# 11/25/2024
# Using chatgpt code copilot

import MAGE
import numpy as np


# --- Generate synthetic data for testing ---
num_genes = 100
num_replicates = 10

# Random gene expression data for X and Y (simulating biological data)
data_x = np.random.rand(num_genes, num_replicates) * 100
data_y = np.random.rand(num_genes, num_replicates) * 100

# --- Run MAGE function ---
print("Running MAGE...")
density_matrix = MAGE.process(data_x, data_y, grid_density=50, output_plots=True)