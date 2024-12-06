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
        contour_info[0], contour_info[2], optimal_cont_lv, 
        grid_width, grid_height, start_x, start_y, 
        gene_mean_x, gene_mean_y, gene_std_x, gene_std_y, num_gene)

    # OS adjustment
    adjusted_score, inliers, outliers = adjust_outlier_scores(
        data_x, data_y, outlier_score, outlier_dist_score, 
        remove_high_low_expr, contour_loops_optimal_info, num_genes)
    
    return adjusted_score, fdr


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
            contours[optimal_lv_idx],
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


def adjust_outlier_scores(data_x, data_y, outlier_score, outlier_dist_score, 
                          remove_high_low_expr, contour_loops_optimal_info, num_genes):
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
        
        containment = contour_loops_optimal_info[-1, 0]
        low_x = mean_data_x_sort[int(np.floor((1 - containment) / 2 * num_genes))]
        low_y = mean_data_y_sort[int(np.floor((1 - containment) / 2 * num_genes))]
        high_x = mean_data_x_sort[int(np.ceil((containment + (1 - containment) / 2) * num_genes))]
        high_y = mean_data_y_sort[int(np.ceil((containment + (1 - containment) / 2) * num_genes))]

        # Remove scores for genes with low/high expression
        for i in range(num_genes):
            mean_x = np.mean(data_x[i, :])
            mean_y = np.mean(data_y[i, :])
            if (mean_x <= low_x and mean_y <= low_y) or (mean_x >= high_x and mean_y >= high_y):
                outlier_score[i] = 1

    # Cap outlier score at 1
    outlier_score[outlier_score > 1] = 1

    # Rank order distance scores < 0
    d_penalty = np.zeros_like(outlier_dist_score)
    negative_d_indices = np.where(outlier_dist_score < 0)[0]
    ranked_indices = np.argsort(np.abs(outlier_dist_score[negative_d_indices]))
    d_penalty[negative_d_indices[ranked_indices]] = np.arange(len(negative_d_indices))
    d_penalty[negative_d_indices] = d_penalty[negative_d_indices] / len(negative_d_indices)

    # Adjust outlier score by standard deviation
    adjusted_outlier_score = 1 - outlier_score - d_penalty
    adjusted_outlier_score[adjusted_outlier_score < 0] = 0

    # Identify genes with top 5% highest outlier scores
    sorted_indices = np.argsort(adjusted_outlier_score)
    in_count = sorted_indices[:int(0.95 * len(adjusted_outlier_score))]
    in_indices = np.zeros_like(adjusted_outlier_score, dtype=bool)
    in_indices[in_count] = True
    out_indices = np.where(~in_indices)[0]

    return adjusted_outlier_score, np.where(in_indices)[0], out_indices


def mage(data_x, data_y, grid_density=50, num_contours=20, output_plots=True, 
         target_containment=0.95, remove_high_low_expr=True, contour_loops_max=10, 
         num_starting_contours=200, monte_carlo_sample_size=1000):
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
    x_coords = start_x + np.arange(grid_density) * grid_width
    y_coords = start_y + np.arange(grid_density) * grid_height

    for k in range(num_genes):
        x_pdf = norm.pdf(x_coords[:, None], gene_mean_x[k], gene_std_x[k])  # Shape: (grid_density, 1)
        y_pdf = norm.pdf(y_coords[None, :], gene_mean_y[k], gene_std_y[k])  # Shape: (1, grid_density)
        density_mat += np.dot(x_pdf, y_pdf)

    # --- Step 3: Determine CER (Characteristic Expression Region) ---
    print("Determining characteristic expression region...")
    density_mat = np.pad(density_mat, pad_width=1, mode='constant', constant_values=0)
    contour_range = np.linspace(np.min(density_mat), np.max(density_mat), num_starting_contours)
    contour_loops_optimal_info = np.zeros((contour_loops_max, 3))

    for c in range(contour_loops_max):
        cs = plt.contour(np.log(density_mat + 1), levels=contour_range)
        contours = cs.allsegs
        levels = cs.levels
        plt.close()  # Close the contour plot if no display is needed

        gene_prob_contained_in_lv = np.zeros(len(levels))

        for lv, contour_set in enumerate(contours):
            if len(contour_set) == 0:  # Skip empty contours
                continue
            path = Path(contour_set[0])  # Take the first region as representative
            monte_carlo_points = np.column_stack([
                np.random.uniform(start_x, start_x + grid_density * grid_width, monte_carlo_sample_size),
                np.random.uniform(start_y, start_y + grid_density * grid_height, monte_carlo_sample_size),
            ])
            gene_prob_contained_in_lv[lv] = np.mean(path.contains_points(monte_carlo_points))

        cont_effectiveness = 1 - np.abs(target_containment - gene_prob_contained_in_lv)
        optimal_lv_idx = np.argmax(cont_effectiveness)
        contour_loops_optimal_info[c] = [gene_prob_contained_in_lv[optimal_lv_idx], 
                                         levels[optimal_lv_idx], 
                                         cont_effectiveness[optimal_lv_idx]]

        if c > 0 and contour_loops_optimal_info[c, 2] <= contour_loops_optimal_info[c - 1, 2]:
            contour_loops_optimal_info = contour_loops_optimal_info[:c]
            break

    optimal_cont_lv = int(contour_loops_optimal_info[-1, 1])

    if output_plots:
        visualize_all_contours(density_mat, contour_range, start_x, start_y, grid_width, grid_height, 
                                        gene_mean_x, gene_mean_y, gene_std_x, gene_std_y)

    # --- Step 4: Assign outlier scores ---
    print("Assigning outlier scores...")
    outlier_score = np.zeros(num_genes)
    outlier_dist_score = np.zeros(num_genes)

    # Vectorized Monte Carlo sampling for outlier scores
    for region in contours[optimal_cont_lv]:
        path = Path(region)
        rand_points = np.column_stack([
            np.repeat(gene_mean_x, monte_carlo_sample_size) + 
            np.repeat(gene_std_x, monte_carlo_sample_size) * np.random.randn(num_genes * monte_carlo_sample_size),
            np.repeat(gene_mean_y, monte_carlo_sample_size) + 
            np.repeat(gene_std_y, monte_carlo_sample_size) * np.random.randn(num_genes * monte_carlo_sample_size),
        ])
        gene_indices = np.repeat(np.arange(num_genes), monte_carlo_sample_size)
        inside_mask = path.contains_points(rand_points)
        inside_gene_counts = np.bincount(gene_indices[inside_mask], minlength=num_genes)
        outlier_score += inside_gene_counts / monte_carlo_sample_size

    # --- Step 5: Adjust outlier scores ---
    print("Adjusting outlier scores...")
    adjusted_score, inliers, outliers = adjust_outlier_scores(
        data_x, data_y, outlier_score, outlier_dist_score, 
        remove_high_low_expr, contour_loops_optimal_info, num_genes
    )

    # --- Step 6: Output Results ---
    if output_plots:
        visualize_outlier_scores(data_x, data_y, adjusted_score, inliers, outliers, 
                             contours, optimal_cont_lv)
        
    return adjusted_score, inliers, outliers


def visualize_all_contours(
    density_mat, contour_range, start_x, start_y, grid_width, grid_height, 
    gene_mean_x, gene_mean_y, gene_std_x=None, gene_std_y=None
):
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

    # Overlay gene positions
    plt.scatter(
        gene_mean_x, gene_mean_y, color="red", s=50, edgecolor="k", alpha=0.7, label="Gene Mean Positions"
    )

    # Optionally, add error bars for gene expression
    if gene_std_x is not None and gene_std_y is not None:
        for i in range(len(gene_mean_x)):
            plt.errorbar(
                gene_mean_x[i], gene_mean_y[i], 
                xerr=gene_std_x[i], yerr=gene_std_y[i], 
                fmt='o', color='red', alpha=0.3, capsize=2
            )

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
    plt.xlabel("Mean Expression (X)")
    plt.ylabel("Mean Expression (Y)")
    plt.title("All Contours with Gene Positions (Pre-CER Selection)")
    plt.colorbar(cs, label="Log-Density Levels")
    plt.legend()
    plt.grid(True)
    plt.show()


def visualize_outlier_scores(data_x, data_y, adjusted_scores, inliers, outliers, 
                             contours, optimal_cont_lv):
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

    # Extract gene means (needed for visualization)
    gene_mean_x = np.mean(data_x, axis=1)
    gene_mean_y = np.mean(data_y, axis=1)

    # Scatter plot of gene means colored by outlier scores
    plt.figure(figsize=(10, 6))
    scatter = plt.scatter(
        gene_mean_x, gene_mean_y, c=adjusted_scores, cmap='coolwarm', s=50, edgecolor='k', alpha=0.7
    )
    plt.colorbar(scatter, label='Outlier Score')
    plt.title('Gene Expression with Outlier Scores')
    plt.xlabel('Mean Expression (X)')
    plt.ylabel('Mean Expression (Y)')
    plt.grid(True)

    # Overlay the contour
    for region in contours[optimal_cont_lv]:
        plt.plot(region[:, 0], region[:, 1], color='black', linestyle='--', label='Optimal Contour')
    plt.legend()
    plt.show()

    # Plot histogram of adjusted outlier scores
    plt.figure(figsize=(8, 6))
    plt.hist(adjusted_scores, bins=30, color='skyblue', edgecolor='black', alpha=0.7)
    plt.axvline(np.percentile(adjusted_scores, 95), color='red', linestyle='--', label='95% Threshold')
    plt.title('Distribution of Adjusted Outlier Scores')
    plt.xlabel('Outlier Score')
    plt.ylabel('Frequency')
    plt.legend()
    plt.show()

    # Scatter plot highlighting inliers and outliers
    plt.figure(figsize=(10, 6))
    plt.scatter(gene_mean_x[inliers], gene_mean_y[inliers], color='green', label='Inliers', alpha=0.7)
    plt.scatter(gene_mean_x[outliers], gene_mean_y[outliers], color='red', label='Outliers', alpha=0.7)
    plt.title('Inliers vs Outliers in Gene Expression')
    plt.xlabel('Mean Expression (X)')
    plt.ylabel('Mean Expression (Y)')
    plt.legend()
    plt.grid(True)
    plt.show()


# Example usage
if __name__ == "__main__":
    # Test data (replace with real data)
    data_x = np.random.rand(100)
    data_y = np.random.rand(100)
    process(data_x, data_y, grid_density=100, num_contours=30)
