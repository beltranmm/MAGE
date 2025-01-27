

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.path import Path
import pickle


# Density Matrix
try:
    with open("dm.pkl", "rb") as f:
        data = pickle.load(f)

    density_mat = data[0]
    units = data[1]

    plt.figure(figsize=(8, 6))
    plt.imshow(density_mat, cmap='hot', interpolation='nearest', origin='lower')
    plt.colorbar(label="PDF")
    plt.title("Density Matrix Heatmap")
    plt.xlabel('Mean Expression (' + units + ')')
    plt.ylabel('Mean Expression (' + units + ')')
    plt.show()

except EOFError:
    print('Density Matrix figure data not found...')



# MAGE OS
try:
    with open("OS.pkl", "rb") as f:
        data = pickle.load(f)

    gene_mean_x = data[0]
    gene_mean_y = data[1]
    adjusted_scores = data[2]
    monte_carlo_points = data[3]
    cer = data[4]
    start_x = data[5]
    start_y = data[6]
    grid_width = data[7]
    grid_height = data[8]
    units = data[9]

    # convert gene means (TPM) points to grid scale
    gene_mean_x_grid = (gene_mean_x-start_x)/grid_width + 2
    gene_mean_y_grid = (gene_mean_y-start_y)/grid_height + 2

    # Scatter plot
    plt.figure(figsize=(10, 6))
    plt.scatter(monte_carlo_points[:,0], monte_carlo_points[:,1],
                 c='g', s=5, edgecolor='g', alpha=0.5, label='MC Sampling')
    scatter = plt.scatter(
        gene_mean_x_grid, gene_mean_y_grid, c=adjusted_scores, cmap='magma', vmin=0, vmax=1,
          s=50, edgecolor='k', alpha=0.7, label='Gene Mean')
    plt.colorbar(scatter, label='Outlier Score')

    # Overlay the CER
    for region in cer:
        plt.plot(region[:, 0], region[:, 1], color='b', linestyle='--', label='CER')

    plt.xlabel('Mean Expression (' + units + ')')
    plt.ylabel('Mean Expression (' + units + ')')
    locs = plt.xticks()[0]
    labels = np.round((locs-2)*grid_width + start_x)
    plt.xticks(locs,labels)
    locs = plt.yticks()[0]
    labels = np.round((locs-2)*grid_width + start_y)
    plt.yticks(locs,labels)
    plt.legend()
    for pos in ['right', 'top']:
        plt.gca().spines[pos].set_visible(False)
    plt.tick_params(direction='in')
    plt.show()

except:
    print('OS figure data not found...')



try:
    with open("fdr.pkl", "rb") as f:
        data = pickle.load(f)

    plt.figure(figsize=(8, 6))
    plt.scatter(data[:, 0], data[:, 1], color='blue')
    plt.title("FDR at Each Outlier Score Cutoff")
    plt.xlabel("Outlier Score Cutoff")
    plt.ylabel("FDR")
    plt.grid(True)
    plt.show()

except:
    print('FDR figure data not found...')