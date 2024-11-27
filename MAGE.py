# Matthew Beltran
# 11/25/2024

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.path import Path
from scipy.stats import norm
from scipy.spatial.distance import euclidean

def process(data_x, data_y, grid_density=50, num_contours=20, output_plots=False, 
         target_containment=0.95, remove_high_low_expr=True, contour_loops_max=10, 
         num_starting_contours=200):
    """
    Parameters:
    - data_x: np.ndarray
        Input array representing the x-dimension of the data.
    - data_y: np.ndarray
        Input array representing the y-dimension of the data.
    - grid_density: int, optional
        Density of the grid used in analysis (default=50).
    - num_contours: int, optional
        Number of contours for computation (default=20).
    - output_plots: bool, optional
        Flag to enable or disable output plots (default=False).
    - target_containment: float, optional
        Target data containment level (default=0.95).
    - remove_high_low_expr: bool, optional
        Flag to remove high and low expression data (default=True).
    - contour_loops_max: int, optional
        Maximum iterations of contour loops (default=10).
    - num_starting_contours: int, optional
        Initial number of contours to evaluate (default=200).

    Returns:
    - outlier_score: np.ndarray
        Computed outlier scores.
    - fdr: np.ndarray
        False discovery rate values.
    """

    # Initialization and checks (placeholder for further computation steps)
    assert len(data_x) == len(data_y), "data_x and data_y must have the same length."
    
    # Replace with MAGE computation logic
    outlier_score = np.zeros(len(data_x))  # Placeholder
    fdr = np.zeros(len(data_x))  # Placeholder



    #calculate gene mean and STD across replicate groups
    
    num_gene = data_y.shape[0]

    # Initialize arrays to store statistics
    gene_mean_x = np.zeros(num_gene)
    gene_mean_y = np.zeros(num_gene)
    gene_std_x = np.zeros(num_gene)
    gene_std_y = np.zeros(num_gene)

    # Calculate mean and standard deviation for each gene
    for k in range(num_gene):
        gene_mean_x[k] = np.mean(data_x[k, :])
        gene_mean_y[k] = np.mean(data_y[k, :])
        gene_std_x[k] = np.std(data_x[k, :]) + 0.0001
        gene_std_y[k] = np.std(data_y[k, :]) + 0.0001



    # fit grid
    print("Fitting a grid")

    if output_plots:
        # Display points
        plt.figure(figsize=(10, 5))
        plt.scatter(gene_mean_x, gene_mean_y, s=25, c='r', marker='.')
        plt.axis('equal')
        plt.title('Normalized Expression')
        plt.xlabel('TPM')
        plt.ylabel('TPM')
        plt.show()

    # Calculate grid parameters
    grid_height = (np.ptp(gene_mean_y) + 4 * np.mean(np.abs(gene_std_y))) / grid_density
    grid_width = (np.ptp(gene_mean_x) + 4 * np.mean(np.abs(gene_std_x))) / grid_density

    density_mat = np.zeros((grid_density, grid_density))  # Initialize density matrix

    start_x = np.min(gene_mean_x) - 2 * np.mean(np.abs(gene_std_x)) + 0.00001
    start_y = np.min(gene_mean_y) - 2 * np.mean(np.abs(gene_std_y)) + 0.00001



    print("Summing gene probabilities")

    # Create grid coordinates
    x_coords = start_x + np.arange(grid_density) * grid_width
    y_coords = start_y + np.arange(grid_density) * grid_height

    # Initialize density matrix
    density_mat = np.zeros((grid_density, grid_density))

    # Compute probability contributions for all genes using broadcasting
    for k in range(len(gene_mean_x)):
        # Gaussian PDFs for the current gene
        x_pdf = norm.pdf(x_coords[:, None], gene_mean_x[k], gene_std_x[k])  # Shape: (grid_density, 1)
        y_pdf = norm.pdf(y_coords[None, :], gene_mean_y[k], gene_std_y[k])  # Shape: (1, grid_density)

        # Accumulate probabilities into the density matrix
        density_mat += np.dot(x_pdf, y_pdf)



    if output_plots:
        print("Plotting density matrix")

        # Display density matrix as a heatmap
        plt.figure(figsize=(8, 8))
        plt.imshow(np.flipud(density_mat), aspect='equal', cmap='viridis')  # Flip the matrix vertically
        plt.colorbar(label='Density')
        plt.title('Density Matrix')
        plt.axis('square')
        plt.xticks([])
        plt.yticks([])
        plt.show()

    # Determine characteristic expression region
    contour_info = select_cer(density_mat, gene_mean_x,
                              gene_mean_y, gene_std_x,
                              gene_std_y, grid_width,
                              grid_height, start_x, start_y,
                              grid_density, output_plots=True)
    print("Optimal Contour Info:")
    print(contour_info)

    # Determine outlier scores
    outlier_scores, outlier_dist_scores = assign_outlier_scores(
        density_cont_points, levels, optimal_cont_lv, 
        grid_width, grid_height, start_x, start_y, 
        gene_mean_x, gene_mean_y, gene_std_x, gene_std_y, num_gene)


    return outlier_score, fdr

def select_cer(density_mat, gene_mean_x, gene_mean_y, gene_std_x, gene_std_y, 
               grid_width, grid_height, start_x, start_y, grid_density, 
               num_starting_contours=200, num_contours=20, contour_loops_max=10, 
               target_containment=0.95, output_plots=True):
    """
    Determines the characteristic expression region (CER) by analyzing density contours.

    Parameters:
    - density_mat: np.ndarray
        The density matrix computed from gene probabilities.
    - gene_mean_x, gene_mean_y: np.ndarray
        Mean expression values for each gene.
    - gene_std_x, gene_std_y: np.ndarray
        Standard deviations of expression values for each gene.
    - grid_width, grid_height: float
        Dimensions of each grid cell.
    - start_x, start_y: float
        Starting coordinates of the grid.
    - grid_density: int
        Number of grid cells along each dimension.
    - num_starting_contours: int, optional
        Initial number of contours to evaluate (default=200).
    - num_contours: int, optional
        Number of contour levels (default=20).
    - contour_loops_max: int, optional
        Maximum number of contour refinement loops (default=10).
    - target_containment: float, optional
        Desired containment probability (default=0.95).
    - output_plots: bool, optional
        Flag to enable or disable visualizations (default=True).

    Returns:
    - contour_loops_optimal_info: np.ndarray
        Optimal contour information (containment, size, effectiveness).
    """

    print("Determining characteristic expression region (CER)")
    density_mat = np.pad(density_mat, pad_width=1, mode='constant', constant_values=0)
    contour_loops_optimal_info = np.zeros((contour_loops_max, 3))
    contour_range = np.linspace(np.min(density_mat), np.max(density_mat), num_starting_contours)

    for c in range(contour_loops_max):
        print(f"Contour loop {c + 1}/{contour_loops_max}")

        # Plot contour
        fig, ax = plt.subplots() if output_plots else (None, None)
        cs = ax.contour(np.log(density_mat + 1), levels=contour_range) if output_plots else plt.contour(np.log(density_mat + 1), levels=contour_range)

        # Extract contour levels and points
        contours = cs.allsegs
        levels = cs.levels
        plt.close(fig)  # Close the plot if only computing data

        gene_prob_contained_in_lv = np.zeros(len(levels))
        
        for lv, contour_set in enumerate(contours):
            if len(contour_set) == 0:  # Skip empty contours
                continue

            gene_is_contained = np.zeros(len(gene_mean_x))
            
            for region in contour_set:
                # Monte Carlo sampling from CPDF
                path = Path(region)
                rand_gene_points = np.column_stack([
                    gene_mean_x + gene_std_x * np.random.randn(len(gene_mean_x)),
                    gene_mean_y + gene_std_y * np.random.randn(len(gene_mean_y))
                ])
                gene_is_contained += path.contains_points(rand_gene_points)

            # Calculate percentage of gene probabilities contained
            gene_prob_contained_in_lv[lv] = np.sum(gene_is_contained) / len(gene_is_contained)

        # Determine contour effectiveness
        cont_effectiveness = 1 - np.abs(target_containment - gene_prob_contained_in_lv)
        optimal_lv_idx = np.argmax(cont_effectiveness)

        # Save optimal contour info
        contour_loops_optimal_info[c] = [
            gene_prob_contained_in_lv[optimal_lv_idx],
            levels[optimal_lv_idx],
            cont_effectiveness[optimal_lv_idx],
        ]

        # Stop early if effectiveness decreases
        if c > 0 and contour_loops_optimal_info[c, 2] <= contour_loops_optimal_info[c - 1, 2]:
            contour_loops_optimal_info = contour_loops_optimal_info[:c]
            break

        # Refine contour range
        optimal_level = levels[optimal_lv_idx]
        lower_bound = max(0, optimal_lv_idx - num_contours // 4)
        upper_bound = min(len(levels) - 1, optimal_lv_idx + num_contours // 4)
        contour_range = levels[lower_bound:upper_bound + 1]

    return contour_loops_optimal_info


def assign_outlier_scores(density_cont_points, levels, optimal_cont_lv, 
                          grid_width, grid_height, start_x, start_y, 
                          gene_mean_x, gene_mean_y, gene_std_x, gene_std_y, 
                          num_gene, monte_carlo_sample_size=1000, output_plots=True):
    """
    Assigns outlier scores and calculates distances for each gene based on contours.

    Parameters:
    - density_cont_points: np.ndarray
        Points defining the density contour lines.
    - levels: np.ndarray
        Contour level values.
    - optimal_cont_lv: int
        Index of the optimal contour level.
    - grid_width, grid_height: float
        Grid cell dimensions.
    - start_x, start_y: float
        Starting coordinates for the grid.
    - gene_mean_x, gene_mean_y: np.ndarray
        Mean gene expression values for x and y dimensions.
    - gene_std_x, gene_std_y: np.ndarray
        Standard deviations of gene expression values.
    - num_gene: int
        Number of genes in the dataset.
    - monte_carlo_sample_size: int, optional
        Number of Monte Carlo samples per gene (default=1000).
    - output_plots: bool, optional
        Flag to enable or disable plotting (default=True).

    Returns:
    - outlier_scores: np.ndarray
        Array of outlier scores for each gene.
    - outlier_dist_scores: np.ndarray
        Array of outlier distance scores for each gene.
    """
    print("Assigning outlier scores")

    # Find points corresponding to the optimal contour level
    confidence_cont_inds = np.where(density_cont_points[0, :] == levels[optimal_cont_lv])[0]
    num_vertices = density_cont_points[1, confidence_cont_inds]

    # Initialize scores
    outlier_scores = np.zeros(num_gene)
    outlier_dist_scores = np.zeros(num_gene)
    first_point_calc = np.ones(num_gene, dtype=bool)

    for region_idx, num_vert in enumerate(num_vertices):
        # Extract contour vertices for the region
        start_idx = confidence_cont_inds[region_idx]
        confidence_cont = np.column_stack((
            density_cont_points[0, start_idx + 1:start_idx + 1 + num_vert],
            density_cont_points[1, start_idx + 1:start_idx + 1 + num_vert]
        ))

        # Scale contour back to the original coordinates
        confidence_cont[:, 0] = (confidence_cont[:, 0] - 1) * grid_width + start_x - grid_width
        confidence_cont[:, 1] = (confidence_cont[:, 1] - 1) * grid_height + start_y - grid_height

        # Monte Carlo sampling for each gene
        for k in range(num_gene):
            # Generate random points around the gene's mean
            rand_points = np.column_stack([
                gene_mean_x[k] + gene_std_x[k] * np.random.randn(monte_carlo_sample_size),
                gene_mean_y[k] + gene_std_y[k] * np.random.randn(monte_carlo_sample_size),
            ])

            # Calculate outlier score
            path = Path(confidence_cont)
            inside = path.contains_points(rand_points)
            outlier_scores[k] += np.sum(inside) / monte_carlo_sample_size

            # Calculate shortest distance from mean point to the contour
            for i in range(len(confidence_cont) - 1):
                segment_start = confidence_cont[i]
                segment_end = confidence_cont[i + 1]

                # Vector math to calculate distance from a point to a line segment
                seg_vec = segment_end - segment_start
                point_vec = np.array([gene_mean_x[k], gene_mean_y[k]]) - segment_start
                proj = np.clip(np.dot(point_vec, seg_vec) / np.dot(seg_vec, seg_vec), 0, 1)
                closest_point = segment_start + proj * seg_vec
                dist_to_seg = euclidean([gene_mean_x[k], gene_mean_y[k]], closest_point)

                # Update shortest distance
                if first_point_calc[k] or dist_to_seg < abs(outlier_dist_scores[k]):
                    first_point_calc[k] = False
                    if path.contains_point([gene_mean_x[k], gene_mean_y[k]]):
                        outlier_dist_scores[k] = -dist_to_seg  # Negative if inside
                    else:
                        outlier_dist_scores[k] = dist_to_seg  # Positive if outside

            # Plot a random gene for visualization
            if output_plots and k == np.random.randint(num_gene):
                plt.scatter(rand_points[:, 0], rand_points[:, 1], s=10, color='red', label='MC Points')
                plt.scatter(gene_mean_x[k], gene_mean_y[k], color='blue', s=50, label='Gene Mean')
                plt.plot(confidence_cont[:, 0], confidence_cont[:, 1], color='black', label='Contour')
                plt.title("Scoring Single Gene")
                plt.xlabel("TPM")
                plt.ylabel("TPM")
                plt.legend()
                plt.show()

    return outlier_scores, outlier_dist_scores

# Example usage
if __name__ == "__main__":
    # Test data (replace with real data)
    data_x = np.random.rand(100)
    data_y = np.random.rand(100)
    process(data_x, data_y, grid_density=100, num_contours=30)
