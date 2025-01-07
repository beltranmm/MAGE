# Matthew Beltran
# 11/25/2024

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.path import Path
from scipy.stats import norm
from scipy.spatial.distance import euclidean

def mage(data_x, data_y, grid_density=50, num_contours=20, output_plots=True, 
         target_containment=0.95, remove_high_low_expr=True, contour_loops_max=10, 
         num_starting_contours=200, monte_carlo_sample_size=1000, output_diags=False, units = 'TPM'):
    """
    MAGE (Minimal Area Gaussian Estimator) workflow to compute outlier scores for genes.

    Parameters:
    - data_x, data_y: np.ndarray
        Gene expression data arrays (genes x replicates).
    - grid_density: int, optional
        Density of the grid (default=50).
    - num_contours: int, optional
        Number of contours for computation (default=20).
    - output_plots: bool, optional
        Flag to enable or disable output plots (default=True).
    - target_containment: float, optional
        Desired containment probability (default=0.95).
    - remove_high_low_expr: bool, optional
        Remove high/low expression genes (default=True).
    - contour_loops_max: int, optional
        Maximum iterations for contour refinement (default=10).
    - num_starting_contours: int, optional
        Number of starting contour levels (default=200).
    - monte_carlo_sample_size: int, optional
        Number of Monte Carlo samples per gene (default=1000).

    Returns:
    - adjusted_outlier_score: np.ndarray
        Final adjusted outlier scores.
    - in_indices: np.ndarray
        Indices of genes classified as "inliers."
    - out_indices: np.ndarray
        Indices of genes classified as "outliers."
    """
    assert data_x.shape == data_y.shape, "data_x and data_y must have the same shape"
    assert grid_density > 0, "Grid density must be a positive integer"
    assert 0 < target_containment < 1, "Target containment must be between 0 and 1"


    num_genes = data_x.shape[0]

    # --- Step 1: Calculate gene statistics ---
    print("Calculating gene statistics...")
    gene_mean_x = np.mean(data_x, axis=1)
    gene_mean_y = np.mean(data_y, axis=1)
    gene_std_x = np.std(data_x, axis=1) + 0.0001
    gene_std_y = np.std(data_y, axis=1) + 0.0001

    # --- Step 2: Fit the grid ---
    print("Fitting the grid...")
    grid_height = (np.ptp(gene_mean_y) + 4 * np.mean(np.abs(gene_std_y))) / grid_density
    grid_width = (np.ptp(gene_mean_x) + 4 * np.mean(np.abs(gene_std_x))) / grid_density
    density_mat = np.zeros((grid_density, grid_density))
    start_x = np.min(gene_mean_x) - 2 * np.mean(np.abs(gene_std_x)) + 0.00001
    start_y = np.min(gene_mean_y) - 2 * np.mean(np.abs(gene_std_y)) + 0.00001

    # Monte Carlo sampling to fill the density matrix
    print("Summing gene probabilities...")
    x_coords = np.linspace(start_x, start_x + grid_density * grid_width, grid_density)
    y_coords = np.linspace(start_y, start_y + grid_density * grid_height, grid_density)

    for k in range(num_genes):
        x_pdf = norm.pdf(x_coords[:, None], gene_mean_x[k], gene_std_x[k])  # Shape: (grid_density, 1)
        y_pdf = norm.pdf(y_coords[None, :], gene_mean_y[k], gene_std_y[k])  # Shape: (1, grid_density)
        density_mat += np.dot(x_pdf, y_pdf)

    if output_diags:
        plt.figure(figsize=(8, 6))
        plt.imshow(density_mat, cmap='hot', interpolation='nearest', origin='lower')
        plt.colorbar(label="PDF")
        plt.title("Density Matrix Heatmap")
        plt.show()

    # --- Step 3: Determine CER (Characteristic Expression Region) ---
    print("Determining characteristic expression region...")
    density_mat = np.pad(density_mat, pad_width=1, mode='constant', constant_values=0)
    contour_range = np.linspace(np.min(density_mat), np.max(density_mat), num_starting_contours)
    cer_effectivness = 0

    if output_diags:
        cer_plot_fig = plt.figure(figsize=(10, 6))
        cer_plot = cer_plot_fig.subplots()

    # monte carlo sampling from each gene
    monte_carlo_points = np.zeros((num_genes*monte_carlo_sample_size,2))
    for i in range(num_genes):
        ind = monte_carlo_sample_size*(i-1)
        monte_carlo_points[range(ind,ind+monte_carlo_sample_size),:] = np.column_stack([
                    np.random.normal(gene_mean_x[i], gene_std_x[i], monte_carlo_sample_size),
                    np.random.normal(gene_mean_y[i], gene_std_y[i], monte_carlo_sample_size),])
    
    # convert MC points (TPM) points to grid scale
    monte_carlo_points[:,0] = (monte_carlo_points[:,0]-start_x)/grid_width + 2
    monte_carlo_points[:,1] = (monte_carlo_points[:,1]-start_y)/grid_height + 2
    
    # Find contour that contains specified number of MC points
    for c in range(contour_loops_max):
        contour_fig = plt.figure()
        contour_ax = contour_fig.subplots()
        cs = contour_ax.contour(np.log(density_mat + 1), levels=contour_range)
        contours = cs.allsegs
        levels = cs.levels
        plt.close(contour_fig)  # Close the contour plot if no display is needed

        gene_prob_contained_in_lv = np.zeros((len(levels),monte_carlo_sample_size*num_genes))

        for lv, contour_set in enumerate(contours):
            if len(contour_set) == 0:  # Skip empty contours
                continue
            for s in range(len(contour_set)):
                path = Path(contour_set[s])
                #gene_prob_contained_in_lv[lv,:]+= np.sum(path.contains_points(monte_carlo_points))
                region_contains = path.contains_points(monte_carlo_points)
                gene_prob_contained_in_lv[lv,region_contains] = 1
                #for i in range(monte_carlo_sample_size):
                    #gene_prob_contained_in_lv[lv,i] = any((gene_prob_contained_in_lv[lv,i],region_contains[i]))

        gene_prob_contained_in_lv = np.sum(gene_prob_contained_in_lv,1)/(monte_carlo_sample_size*num_genes)
        cont_effectiveness = 1 - np.abs(target_containment - gene_prob_contained_in_lv)
        optimal_lv_idx = np.argmax(cont_effectiveness)
        
        # exit loop if no improvement
        if c > 0 and cont_effectiveness[optimal_lv_idx] <= cer_effectivness:
            break
        else:
            cer = contours[optimal_lv_idx]
            cer_effectivness = cont_effectiveness[optimal_lv_idx]
            if output_diags:
                for region in cer:
                    cer_plot.plot(region[:, 0], region[:, 1], color='black')

            # select new range of contour levels
            if optimal_lv_idx > 1 and optimal_lv_idx < num_starting_contours:
                contour_range = np.linspace(contour_range[optimal_lv_idx - 1], contour_range[optimal_lv_idx + 1], num_starting_contours)
            else:
                break

    if output_diags:
        cer_plot.set_xlabel('Mean Expression (' + units + ')')
        cer_plot.set_ylabel('Mean Expression (' + units + ')')
        cer_plot.grid(True)
        cer_plot.legend()
        cer_plot_fig.show()
        visualize_all_contours(density_mat, contour_range, start_x, start_y, grid_width, grid_height, 
                                gene_mean_x, gene_mean_y, gene_std_x, gene_std_y, units)

    # --- Step 4: Assign outlier scores ---
    print("Assigning outlier scores...")
    outlier_score = np.zeros(num_genes)
    outlier_dist_score = np.zeros(num_genes)


    # Vectorized Monte Carlo sampling for outlier scores
    outlier_score = np.zeros(num_genes)
    for i in range(num_genes):
        for region in cer:
            if len(region) == 0:
                continue

            path = Path(region)
            ind = monte_carlo_sample_size*(i-1)
            outlier_score[i]+= np.sum(path.contains_points(monte_carlo_points[range(ind,ind+monte_carlo_sample_size),:]))

        
    outlier_score = 1 - outlier_score / monte_carlo_sample_size

    # --- Step 5: Adjust outlier scores ---
    print("Adjusting outlier scores...")
    adjusted_score, inliers, outliers = adjust_outlier_scores(data_x, data_y, outlier_score, outlier_dist_score, 
                                                            remove_high_low_expr, target_containment, num_genes)

    # --- Step 6: Output Results ---
    if output_plots:
        visualize_MC_sampling(gene_mean_x, gene_mean_y, adjusted_score, monte_carlo_points,
                              cer, start_x, start_y, grid_width, grid_height, units)
        
    return adjusted_score


def adjust_outlier_scores(data_x, data_y, outlier_score, outlier_dist_score, 
                          remove_high_low_expr, target_containment, num_genes):
    """
    Adjusts the outlier scores by removing high/low expression outliers and ranking distance scores.

    Parameters:
    - data_x, data_y: np.ndarray
        Input gene expression data arrays (genes x replicates).
    - outlier_score: np.ndarray
        Initial outlier scores.
    - outlier_dist_score: np.ndarray
        Initial outlier distance scores.
    - remove_high_low_expr: bool
        Flag to enable/disable high/low expression removal.
    - contour_loops_optimal_info: np.ndarray
        Optimal contour information from the MAGE process.
    - num_genes: int
        Number of genes.

    Returns:
    - adjusted_outlier_score: np.ndarray
        Adjusted outlier scores.
    - in_indices: np.ndarray
        Indices of genes classified as "inlier."
    - out_indices: np.ndarray
        Indices of genes classified as "outlier."
    """
    print("Adjusting outlier scores")
    
    if remove_high_low_expr:
        print("Removing low/high expression outlier scores")

        # Sort the mean values of data_x and data_y
        mean_data_x_sort = np.sort(np.mean(data_x, axis=1))
        mean_data_y_sort = np.sort(np.mean(data_y, axis=1))
        
        low_x = mean_data_x_sort[int(np.floor((1 - target_containment) / 2 * num_genes))]
        low_y = mean_data_y_sort[int(np.floor((1 - target_containment) / 2 * num_genes))]
        high_x = mean_data_x_sort[int(np.ceil((target_containment + (1 - target_containment) / 2) * num_genes))]
        high_y = mean_data_y_sort[int(np.ceil((target_containment + (1 - target_containment) / 2) * num_genes))]

        # Remove scores for genes with low/high expression
        for i in range(num_genes):
            mean_x = np.mean(data_x[i, :])
            mean_y = np.mean(data_y[i, :])
            if (mean_x <= low_x and mean_y <= low_y) or (mean_x >= high_x and mean_y >= high_y):
                outlier_score[i] = 0

    # Cap outlier score at 1
    outlier_score[outlier_score > 1] = 1

    # Rank order distance scores < 0
    d_penalty = np.zeros_like(outlier_dist_score)
    negative_d_indices = np.where(outlier_dist_score < 0)[0]
    ranked_indices = np.argsort(np.abs(outlier_dist_score[negative_d_indices]))
    d_penalty[negative_d_indices[ranked_indices]] = np.arange(len(negative_d_indices))
    d_penalty[negative_d_indices] = d_penalty[negative_d_indices] / len(negative_d_indices)

    # Adjust outlier score by standard deviation
    adjusted_outlier_score = outlier_score + d_penalty
    adjusted_outlier_score[adjusted_outlier_score < 0] = 0
    adjusted_outlier_score[adjusted_outlier_score > 1] = 1

    # Identify genes with top 5% highest outlier scores
    sorted_indices = np.argsort(adjusted_outlier_score)
    in_count = sorted_indices[:int(0.95 * len(adjusted_outlier_score))]
    in_indices = np.zeros_like(adjusted_outlier_score, dtype=bool)
    in_indices[in_count] = True
    out_indices = np.where(~in_indices)[0]

    return adjusted_outlier_score, np.where(in_indices)[0], out_indices

def FDR(data_x, data_y, outlier_score, grid_density=50, num_contours=20, output_plots=True, 
         target_containment=0.95, remove_high_low_expr=True, contour_loops_max=10, 
         num_starting_contours=200, monte_carlo_sample_size=1000, output_diags=False):
    # permutate data
    data_x, data_y = permute_and_filter(data_x, data_y, outlier_score, 10)

    # call MAGE w/ permutated data
    OS_perm = mage(data_x, data_y, grid_density, num_contours, output_plots,
                    target_containment, remove_high_low_expr, contour_loops_max, 
                    num_starting_contours, monte_carlo_sample_size, output_diags)
    
    # calculate FDR
    FDR = calculate_fdr(outlier_score, OS_perm, 100, output_plots)

def permute_and_filter(dataX, dataY, OutlierScore, num_permutations=10):
    """
    Permute columns of data, and filter out genes with the highest 10% OutlierScore.

    Parameters:
        dataX (ndarray): First data matrix.
        dataY (ndarray): Second data matrix.
        OutlierScore (ndarray): Array containing outlier scores for genes.
        num_permutations (int): Number of times to permute columns. Default is 10.

    Returns:
        dataX_filtered (ndarray): Filtered dataX after removing top 10% outlier genes.
        dataY_filtered (ndarray): Filtered dataY after removing top 10% outlier genes.
    """
    # Step 1: Concatenate dataX and dataY horizontally
    dataTotal = np.hstack((dataX, dataY))

    # Step 2: Permute columns randomly, 'num_permutations' times
    for _ in range(num_permutations):
        permuted_indices = np.random.permutation(dataTotal.shape[1])
        dataTotal = dataTotal[:, permuted_indices]

    # Step 3: Sort indices by OutlierScore (ascending) and keep 90% of genes
    sorted_indices = np.argsort(OutlierScore)  # Sort in ascending order
    keep_indices = sorted_indices[:round(0.9 * len(sorted_indices))]  # Keep bottom 90%

    # Step 4: Filter rows (genes) based on indices
    dataTotal_filtered = dataTotal[keep_indices, :]

    # Step 5: Split back into dataX and dataY
    dataX_filtered = dataTotal_filtered[:, :dataX.shape[1]]
    dataY_filtered = dataTotal_filtered[:, dataX.shape[1]:]

    return dataX_filtered, dataY_filtered

def calculate_fdr(OutlierScore, OutlierScore_perm, num_steps=100, output_plot=False):
    """
    Calculate the False Discovery Rate (FDR) at each Outlier Score (OS).

    Parameters:
        OutlierScore (ndarray): Array of observed outlier scores.
        OutlierScore_perm (ndarray): Array of permuted outlier scores.
        num_steps (int): Number of OS thresholds to evaluate. Default is 100.

    Returns:
        OS_FDR (ndarray): 2D array with OS thresholds and corresponding FDR values.
    """
    # Initialize OS_FDR array
    OS_FDR = np.zeros((num_steps, 2))

    # Initial outlier score cutoff
    outlier_score_cutoff = 0

    # Step size for outlier score cutoff
    step_size = 1 / num_steps

    # Calculate FDR for each threshold
    for i in range(num_steps):
        num_above = np.sum(OutlierScore >= outlier_score_cutoff)
        num_above_perm = np.sum(OutlierScore_perm >= outlier_score_cutoff)

        OS_FDR[i, 0] = outlier_score_cutoff
        OS_FDR[i, 1] = num_above_perm / (num_above + num_above_perm) if (num_above + num_above_perm) > 0 else 0

        outlier_score_cutoff += step_size

    # Replace NaN values in FDR column with 0
    OS_FDR[np.isnan(OS_FDR[:, 1]), 1] = 0

    # Plot FDR scatter plot
    if output_plot:
        plt.figure(figsize=(8, 6))
        plt.scatter(OS_FDR[:, 0], OS_FDR[:, 1], color='blue', label='FDR')
        plt.title("FDR at Each Outlier Score Cutoff")
        plt.xlabel("Outlier Score Cutoff")
        plt.ylabel("FDR")
        plt.grid(True)
        plt.legend()
        plt.show()

    # Assign FDR to genes
    gene_FDR = np.zeros((len(OutlierScore), 1))
    for i in range(len(OutlierScore)):
        ind = np.argmin((abs(OS_FDR[:,0]-OutlierScore[i])))
        gene_FDR[i] = OS_FDR[ind,1]


    return gene_FDR

def visualize_all_contours(density_mat, contour_range, start_x, start_y, grid_width, grid_height, 
    gene_mean_x, gene_mean_y, gene_std_x=None, gene_std_y=None, units = 'TPM'):
    """
    Visualizes all contours from the density matrix and overlays gene positions.

    Parameters:
    - density_mat: np.ndarray
        The density matrix representing gene expression probabilities.
    - contour_range: np.ndarray
        Range of contour levels used to generate the contours.
    - start_x, start_y: float
        Starting coordinates of the grid.
    - grid_width, grid_height: float
        Dimensions of each grid cell.
    - gene_mean_x, gene_mean_y: np.ndarray
        Mean expression values of the genes.
    - gene_std_x, gene_std_y: np.ndarray, optional
        Standard deviations of the genes, for optional error bars.
    """
    print("Visualizing all contours with gene positions...")

    # Create the contour plot
    plt.figure(figsize=(10, 8))
    cs = plt.contour(
        np.log(density_mat + 1), levels=contour_range, cmap="viridis"
    )

    # convert gene means (TPM) points to grid scale
    gene_mean_x = (gene_mean_x-start_x)/grid_width + 2
    gene_mean_y = (gene_mean_y-start_y)/grid_height + 2
    
    # Overlay gene positions
    plt.scatter(
        gene_mean_x, gene_mean_y, color="red", s=50, edgecolor="k", alpha=0.7, label="Gene Mean Positions"
    )

    # Optionally, add error bars for gene expression
    if len(gene_mean_x) < 200:
        if gene_std_x is not None and gene_std_y is not None:
            for i in range(len(gene_mean_x)):
                plt.errorbar(
                    gene_mean_x[i], gene_mean_y[i], 
                    xerr=gene_std_x[i], yerr=gene_std_y[i], 
                    fmt='o', color='red', alpha=0.3, capsize=2)

    # Rescale axes to match original data scale
    plt.xticks(
        ticks=np.linspace(0, density_mat.shape[1], 5),
        labels=[f"{start_x + i * grid_width:.2f}" for i in range(5)],
    )
    plt.yticks(
        ticks=np.linspace(0, density_mat.shape[0], 5),
        labels=[f"{start_y + i * grid_height:.2f}" for i in range(5)],
    )

    # Plot labels and legend
    plt.xlabel('Mean Expression (' + units + ')')
    plt.ylabel('Mean Expression (' + units + ')')
    plt.title("All Contours with Gene Positions (Pre-CER Selection)")
    plt.colorbar(cs, label="Log-Density Levels")
    plt.legend()
    plt.grid(True)
    plt.show()

def visualize_MC_sampling(gene_mean_x, gene_mean_y, adjusted_scores, monte_carlo_points, cer, start_x, start_y, grid_width, grid_height, units):
    """
    Visualizes the outlier scores and characteristic expression regions.

    Parameters:
    - data_x, data_y: np.ndarray
        Input gene expression data arrays (genes x replicates).
    - adjusted_scores: np.ndarray
        Adjusted outlier scores for visualization.
    - inliers, outliers: np.ndarray
        Indices of inliers and outliers.
    - gene_mean_x, gene_mean_y: np.ndarray
        Mean expression values for each gene.
    - contours: list
        List of contour paths from the MAGE process.
    - optimal_cont_lv: int
        Index of the optimal contour level.
    """

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
        plt.plot(region[:, 0], region[:, 1], color='black', linestyle='--', label='CER')

    plt.xlabel('Mean Expression (' + units + ')')
    plt.ylabel('Mean Expression (' + units + ')')
    plt.grid(True)
    plt.legend()
    plt.show()


# Example usage
if __name__ == "__main__":
    # Test data
    data_x = np.random.rand(100,5)
    data_y = np.random.rand(100,5)
    mage(data_x, data_y)
