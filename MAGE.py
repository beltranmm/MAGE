# Matthew Beltran
# 11/25/2024

import numpy as np
import random
import matplotlib.pyplot as plt
from matplotlib.path import Path
from scipy.stats import norm
from scipy.spatial.distance import euclidean
import math
import csv
import pickle

def mage(data_x, data_y, grid_density=50, num_contours=5, output_plots=True, 
         target_containment=0.95, remove_high_low_expr=True, dist_adj=True, contour_loops_max=20, 
         num_starting_contours=20, monte_carlo_sample_size=1000, nonnegative = True, output_diags=False, units = 'TPM',
         notifications = True, saveFigs = False):
    """
    MAGE (Monte-carlo method for Aberrant Gene Expression) workflow to compute outlier scores for genes.

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
    - outlier_score: np.ndarray
        Final outlier scores.
    """
    assert data_x.shape[0] == data_y.shape[0], "data_x and data_y must have the same shape"
    assert grid_density > 0, "Grid density must be a positive integer"
    assert 0 < target_containment < 1, "Target containment must be between 0 and 1"


    num_genes = data_x.shape[0]

    # --- Step 1: Calculate gene statistics ---
    if notifications: print("Calculating gene statistics...")
    gene_mean_x = np.mean(data_x, axis=1)
    gene_mean_y = np.mean(data_y, axis=1)
    gene_std_x = np.std(data_x, axis=1) + 0.0001
    gene_std_y = np.std(data_y, axis=1) + 0.0001

    # --- Step 2: Fit the grid ---
    if notifications: print("Fitting the grid...")
    grid_height = (np.ptp(gene_mean_y) + 4 * np.mean(np.abs(gene_std_y))) / grid_density
    grid_width = (np.ptp(gene_mean_x) + 4 * np.mean(np.abs(gene_std_x))) / grid_density
    density_mat = np.zeros((grid_density, grid_density))
    start_x = np.min(gene_mean_x) - 2 * np.mean(np.abs(gene_std_x)) + 0.00001
    start_y = np.min(gene_mean_y) - 2 * np.mean(np.abs(gene_std_y)) + 0.00001

    # Gaussian PDF to fill the density matrix
    if notifications: print("Summing gene probabilities...")
    x_coords = np.linspace(start_x, start_x + grid_density * grid_width, grid_density)
    y_coords = np.linspace(start_y, start_y + grid_density * grid_height, grid_density)

    for k in range(num_genes):
        x_pdf = norm.pdf(x_coords[:, None], gene_mean_x[k], gene_std_x[k])  # Shape: (grid_density, 1)
        y_pdf = norm.pdf(y_coords[None, :], gene_mean_y[k], gene_std_y[k])  # Shape: (1, grid_density)
        density_mat += np.dot(x_pdf, y_pdf)

    if nonnegative:
        for i in range(grid_density):
            for j in range(grid_density):
                if x_coords[i] <= 0 or y_coords[j] <= 0:
                    density_mat[i,j] = 0

    if saveFigs:
        with open("dm.pkl", "wb") as f:
            pickle.dump((density_mat, units), f)

    if output_diags:
        plt.figure(figsize=(8, 6))
        plt.imshow(density_mat, cmap='hot', interpolation='nearest', origin='lower')
        plt.colorbar(label="PDF")
        plt.title("Density Matrix Heatmap")
        plt.xlabel('Mean Expression (' + units + ')')
        plt.ylabel('Mean Expression (' + units + ')')
        plt.show()

    # transpose density matrix
    density_mat = density_mat.transpose()

    # --- Step 3: Determine CER (Characteristic Expression Region) ---
    if notifications: print("Determining characteristic expression region...")
    density_mat = np.pad(density_mat, pad_width=1, mode='constant', constant_values=0)
    density_mat = np.log(density_mat + 1)
    contour_range_min = np.min(density_mat)
    contour_range_max = np.max(density_mat)
    contour_range = np.linspace(contour_range_min, contour_range_max, num_starting_contours)
    contour_range_history = contour_range
    cer_effectivness = 0

    if output_diags:
        cer_plot_fig = plt.figure(figsize=(10, 6))
        cer_plot = cer_plot_fig.subplots()

    # get MC points
    monte_carlo_sample_size_CER = round(0.1*monte_carlo_sample_size)
    monte_carlo_points = MC_points(gene_mean_x,gene_mean_y,gene_std_x,gene_std_y,monte_carlo_sample_size_CER,nonnegative)
    # convert MC points (TPM) points to grid scale
    monte_carlo_points[:,0] = (monte_carlo_points[:,0]-start_x)/grid_width + 2
    monte_carlo_points[:,1] = (monte_carlo_points[:,1]-start_y)/grid_height + 2
    
    # Find contour that contains specified number of MC points
    for c in range(contour_loops_max):
        contour_fig = plt.figure()
        contour_ax = contour_fig.subplots()
        contour_range = np.sort(contour_range)
        cs = contour_ax.contour(density_mat, levels=contour_range)
        contours = cs.allsegs
        levels = cs.levels
        plt.close(contour_fig)  # Close the contour plot if no display is needed

        gene_prob_contained_in_lv = np.zeros(len(levels))

        for lv, contour_set in enumerate(contours):
            if len(contour_set) == 0:  # Skip empty contours
                continue
            gene_prob_contained = np.zeros(monte_carlo_sample_size_CER*num_genes)
            for s in range(len(contour_set)):
                path = Path(contour_set[s])
                region_contains = path.contains_points(monte_carlo_points)
                gene_prob_contained[region_contains] = 1
            gene_prob_contained_in_lv[lv] = np.sum(gene_prob_contained)/(monte_carlo_sample_size_CER*num_genes)
                

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
            if optimal_lv_idx > 1:
                contour_range_min = contour_range[optimal_lv_idx - 1]
            if optimal_lv_idx < num_starting_contours:
                contour_range_max = contour_range[optimal_lv_idx + 1]
            
            contour_range = np.linspace(contour_range_min, contour_range_max, num_contours)
            contour_range_history = np.append(contour_range_history,contour_range)
            num_starting_contours = num_contours
    
    if output_diags:
        cer_plot.set_xlabel('Mean Expression (' + units + ')')
        cer_plot.set_ylabel('Mean Expression (' + units + ')')
        cer_plot.grid(True)
        cer_plot.legend()
        cer_plot_fig.show()
        contour_range_history = np.unique(contour_range_history)
        contour_range_history = np.sort(contour_range_history)
        visualize_all_contours(density_mat, contour_range_history, start_x, start_y, grid_width, grid_height, 
                                gene_mean_x, gene_mean_y, gene_std_x, gene_std_y, units)

    # --- Step 4: Assign outlier scores ---
    if notifications: print("Assigning outlier scores...")
    outlier_score = np.zeros(num_genes)
    outlier_dist_score = np.zeros(num_genes)

    # get new MC points
    monte_carlo_points = MC_points(gene_mean_x,gene_mean_y,gene_std_x,gene_std_y,monte_carlo_sample_size,nonnegative)
    # convert MC points and gene means (TPM) points to grid scale
    monte_carlo_points[:,0] = (monte_carlo_points[:,0]-start_x)/grid_width + 2
    monte_carlo_points[:,1] = (monte_carlo_points[:,1]-start_y)/grid_height + 2
    gene_mean_x_grid = (gene_mean_x-start_x)/grid_width + 2
    gene_mean_y_grid = (gene_mean_y-start_y)/grid_height + 2

    # Vectorized Monte Carlo sampling for outlier scores
    for i in range(num_genes):
        for region in cer:
            if len(region) == 0:
                continue

            path = Path(region)
            ind = monte_carlo_sample_size*(i-1)
            outlier_score[i]+= np.sum(path.contains_points(monte_carlo_points[range(ind,ind+monte_carlo_sample_size),:]))
            dist = distance_point_to_path((gene_mean_x_grid[i],gene_mean_y_grid[i]),region)
            if outlier_dist_score[i] == 0 or abs(outlier_dist_score[i]) > dist:
                # check if inside/outside CER (inside is negative)
                if path.contains_point((gene_mean_x_grid[i],gene_mean_y_grid[i])):
                    outlier_dist_score[i] = -1*dist
                else:
                    outlier_dist_score[i] = dist
        
    outlier_score = 1 - outlier_score / monte_carlo_sample_size

    # --- Step 5: Adjust outlier scores ---
    if notifications: print("Adjusting outlier scores...")
    adjusted_score = adjust_outlier_scores(data_x, data_y, outlier_score, outlier_dist_score, 
                                                            remove_high_low_expr, dist_adj,
                                                            target_containment, num_genes,
                                                            notifications=notifications)

    # --- Step 6: Output Results ---
    if output_plots:
        output_figure(gene_mean_x, gene_mean_y, adjusted_score, monte_carlo_points,
                              cer, start_x, start_y, grid_width, grid_height, units)
        
    if saveFigs:
        with open("OS.pkl", "wb") as f:
            pickle.dump((gene_mean_x, gene_mean_y, adjusted_score, monte_carlo_points,
                              cer, start_x, start_y, grid_width, grid_height, units), f)
        
    return adjusted_score

def MC_points(mean_x, mean_y, std_x, std_y, sample_size, nonnegative=True):
    # monte carlo sampling from each gene
    num_gene = len(mean_x)
    points = np.zeros((num_gene*sample_size,2))
    if nonnegative:
        for i in range(num_gene):
            keep_points = np.zeros((sample_size,2))
            neg = True
            while neg:
                indX = np.asarray(keep_points[:,0] <= 0).nonzero()[0]
                indY = np.asarray(keep_points[:,1] <= 0).nonzero()[0]
                if len(indX) == 0 and len(indY) == 0:
                    neg = False
                else:
                    if len(indX) != 0:
                        new_pointsX = np.random.normal(mean_x[i], std_x[i], len(indX))
                        keep_points[indX,0] = new_pointsX

                    if len(indY) != 0:
                        new_pointsY = np.random.normal(mean_y[i], std_y[i], len(indY))
                        keep_points[indY,1] = new_pointsY
            
            ind = sample_size*(i-1)
            points[range(ind,ind+sample_size),:] = keep_points
    else:
        for i in range(num_gene):
            ind = sample_size*(i-1)
            points[range(ind,ind+sample_size),:] = np.column_stack([
                        np.random.normal(mean_x[i], std_x[i], sample_size),
                        np.random.normal(mean_y[i], std_y[i], sample_size),])

    return points

def adjust_outlier_scores(data_x, data_y, outlier_score, outlier_dist_score, 
                          remove_high_low_expr, dist_adj, target_containment, num_genes, notifications = True):
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
    if notifications: print("Adjusting outlier scores")
    
    if remove_high_low_expr:
        if notifications: print("Removing low/high expression outlier scores")

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

    if dist_adj:
        # Rank order distance scores < 0
        d_penalty = np.zeros_like(outlier_dist_score)
        negative_d_indices = np.where(outlier_dist_score < 0)[0]
        ranked_indices = np.argsort(np.abs(outlier_dist_score[negative_d_indices]))
        d_penalty[negative_d_indices[ranked_indices]] = np.arange(len(negative_d_indices))
        d_penalty[negative_d_indices] = d_penalty[negative_d_indices] / len(negative_d_indices)

        # Adjust outlier score by standard deviation
        adjusted_outlier_score = outlier_score - d_penalty
        adjusted_outlier_score[adjusted_outlier_score < 0] = 0
        adjusted_outlier_score[adjusted_outlier_score > 1] = 1
    else:
        adjusted_outlier_score = outlier_score

    return adjusted_outlier_score

def FDR(data_x, data_y, outlier_score, grid_density=50, num_contours=5, output_plots=True, 
         target_containment=0.95, remove_high_low_expr=True, dist_adj=True, contour_loops_max=20, 
         num_starting_contours=20, monte_carlo_sample_size=1000, nonnegative=True, output_diags=False, saveFigs=False):
    
    # permutate data
    data_x, data_y = permute_and_filter(data_x, data_y, outlier_score, 10)

    # call MAGE w/ permutated data
    OS_perm = mage(data_x, data_y, grid_density, num_contours, output_plots,
                    target_containment, remove_high_low_expr, dist_adj, contour_loops_max, 
                    num_starting_contours, monte_carlo_sample_size, nonnegative, output_diags)
    
    # calculate FDR
    FDR = calculate_fdr(outlier_score, OS_perm, 100, output_plots, saveFigs)

    return np.squeeze(FDR)

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

def calculate_fdr(OutlierScore, OutlierScore_perm, num_steps=100, output_plot=False, saveFigs = False):
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

    if saveFigs:
        with open("fdr.pkl", "wb") as f:
            pickle.dump(OS_FDR, f)

    # Assign FDR to genes
    gene_FDR = np.zeros((len(OutlierScore), 1))
    for i in range(len(OutlierScore)):
        ind = np.argmin((abs(OS_FDR[:,0]-OutlierScore[i])))
        gene_FDR[i] = OS_FDR[ind,1]


    return gene_FDR

def distance_point_to_segment(point, segment_start, segment_end):
    """Calculates the distance from a point to a line segment."""
    
    p_x, p_y = point
    s_x1, s_y1 = segment_start
    s_x2, s_y2 = segment_end
    
    segment_vector_x = s_x2 - s_x1
    segment_vector_y = s_y2 - s_y1
    
    point_vector_x = p_x - s_x1
    point_vector_y = p_y - s_y1
    
    segment_length_squared = segment_vector_x**2 + segment_vector_y**2
    
    if segment_length_squared == 0:
      return math.sqrt(point_vector_x**2 + point_vector_y**2)
    
    t = max(0, min(1, (point_vector_x * segment_vector_x + point_vector_y * segment_vector_y) / segment_length_squared))
    
    closest_x = s_x1 + t * segment_vector_x
    closest_y = s_y1 + t * segment_vector_y
    
    return math.sqrt((p_x - closest_x)**2 + (p_y - closest_y)**2)

def distance_point_to_path(point, path):
    """Calculates the shortest distance from a point to a path."""
    
    min_distance = float('inf')
    
    for i in range(len(path) - 1):
        segment_start = path[i]
        segment_end = path[i+1]
        
        distance = distance_point_to_segment(point, segment_start, segment_end)
        min_distance = min(min_distance, distance)
    
    return min_distance

def analyze_depth(dataX, dataY, depths, units = 'TPM', top_cutoff = 0.05, replacement=False, save_data=False):

    OS = np.zeros((dataX.shape[0],len(depths)))

    if save_data:
        dataSave = np.zeros((dataX.shape[0], dataX.shape[1], 2, len(depths)))

    for d in range(len(depths)):
        print("Running depth " + str(d+1) + " of " + str(len(depths)))
        dataX_adj, dataY_adj = adjust_depth(dataX, dataY, depths[d], replacement)

        if save_data:
            dataSave[:,:,0,d] = dataX_adj
            dataSave[:,:,1,d] = dataY_adj

        OS[:,d] = mage(dataX_adj, dataY_adj, output_plots= False, units= units, notifications=False, remove_high_low_expr=False)

    if save_data:
        np.savetxt('depth_data.csv', dataSave, fmt='%s', delimiter=',')

    OS_original = mage(dataX, dataY, output_plots=False, units=units, notifications=False, remove_high_low_expr=False)
    sort_ind = np.argsort(OS_original)
    top_ind_orig = sort_ind[-int(np.ceil(len(OS_original)*top_cutoff)):]

    agreement = np.zeros((len(depths),1))

    for d in range(len(depths)):
        sort_ind = np.argsort(OS[:,d])
        top_ind = sort_ind[-int(np.ceil(OS.shape[0]*top_cutoff)):]
        overlap = np.intersect1d(top_ind,top_ind_orig)
        agreement[d] = len(overlap)/int(np.ceil(len(OS_original)*top_cutoff))

    plt.figure(figsize=(10, 8))
    plt.plot(depths, agreement, color="red")
    plt.scatter(depths, agreement, color="red", s=50)
    plt.xlabel('Sampling Depth %')
    plt.ylabel('% Agreement')
    plt.show()

    return agreement

def adjust_depth(dataX, dataY, depth, replacement=False):
    # convert to integers
    num_gene = dataX.shape[0]
    dataX = dataX.round()
    dataY = dataY.round()


    # est. reads per gene per replicate set

    for rep in range(dataX.shape[1]):
        totalReadX = np.sum(dataX[:,rep])
        totalReadY = np.sum(dataY[:,rep])

        listedReadX = np.zeros((int(totalReadX),1))
        ind = 0
        for i in range(num_gene):
            prevInd = ind
            ind = ind + int(dataX[i,rep])
            listedReadX[prevInd:ind] = i

        listedReadY = np.zeros((int(totalReadY),1))
        ind = 0
        for i in range(num_gene):
            prevInd = ind
            ind = ind + int(dataY[i,rep])
            listedReadY[prevInd:ind] = i

        # Binary sampling of indices
        if replacement or depth > 1:
            sampledReadListX = random.choices(listedReadX.tolist(),k=int(depth*totalReadX))
            sampledReadListY = random.choices(listedReadY.tolist(),k=int(depth*totalReadY))
        else:
            sampledReadListX = random.sample(listedReadX.tolist(), int(depth*totalReadX))
            sampledReadListY = random.sample(listedReadY.tolist(), int(depth*totalReadY))

        # index list
        sampledReadX = [0]*len(sampledReadListX)
        sampledReadY = [0]*len(sampledReadListY)
        for i in range(len(sampledReadListX)):
            sampledReadX[i] = sampledReadListX[i][0]
        for i in range(len(sampledReadListY)):
            sampledReadY[i] = sampledReadListY[i][0]
    


        # calculate new expression with reduced depth
        for i in range(num_gene):
            dataX[i,rep] = sampledReadX.count(i)*(1/depth)
            dataY[i,rep] = sampledReadY.count(i)*(1/depth)
    
    return dataX, dataY

def analyze_samples(dataX, dataY, samples, trials=0, top_cutoff = 0.05, replacement=False, saveData=False):

    OS_original = mage(dataX, dataY, output_plots=False, notifications=False, remove_high_low_expr=False)
    sort_ind = np.argsort(OS_original)
    top_ind_orig = sort_ind[-int(np.ceil(len(OS_original)*top_cutoff)):]

    agreement = np.zeros((len(samples),trials))
    std_dev = np.zeros(len(samples))

    for t in range(trials):
        if trials > 0: print("Running trial " + str(t+1) + " of " + str(trials))

        OS = np.zeros((dataX.shape[0],len(samples)))

        for d in range(len(samples)):
            print("Running sampling " + str(d+1) + " of " + str(len(samples)))
            dataX_adj, dataY_adj = adjust_samples(dataX, dataY, samples[d], replacement)
            OS[:,d] = mage(dataX_adj, dataY_adj, output_plots= False, notifications=False, remove_high_low_expr=False)

    


        for d in range(len(samples)):
            sort_ind = np.argsort(OS[:,d])
            top_ind = sort_ind[-int(np.ceil(OS.shape[0]*top_cutoff)):]
            overlap = np.intersect1d(top_ind,top_ind_orig)
            agreement[d,t] = len(overlap)/int(np.ceil(len(OS_original)*top_cutoff))

    for d in range(len(samples)):
        std_dev[d] = np.std(agreement[d])

    agreement = np.mean(agreement,1)

    if saveData:
        with open("sampleDepth.pkl", "wb") as f:
            pickle.dump((samples, agreement, std_dev), f)

    plt.figure(figsize=(10, 8))
    plt.plot(samples, agreement, color="red")
    plt.scatter(samples, agreement, color="red", s=50)
    plt.errorbar(samples, agreement, yerr=std_dev, fmt='none', ecolor='r')
    plt.ylim(0,1)
    plt.xlabel('Sampling %')
    plt.ylabel('% Agreement')
    plt.show()

    return agreement

def adjust_samples(dataX, dataY, samples, replacement=False):
    
    num_sampleX = dataX.shape[1]
    num_sampleY = dataY.shape[1]

    # Binary sampling
    if replacement or samples > 1:
        Xind = random.choices(range(num_sampleX),k=round(int(samples*num_sampleX)))
        Yind = random.choices(range(num_sampleY),k=round(int(samples*num_sampleY)))
    else:
        Xind = random.sample(range(num_sampleX), round(int(samples*num_sampleX)))
        Yind = random.sample(range(num_sampleY), round(int(samples*num_sampleY)))

    # index list
    dataX_adj = np.zeros((dataX.shape[0],len(Xind)))
    dataY_adj = np.zeros((dataY.shape[0],len(Yind)))

    for i,ind in enumerate(Xind):
        dataX_adj[:,i] = dataX[:,ind]
    for i,ind in enumerate(Yind):
        dataY_adj[:,i] = dataY[:,ind]
    
    return dataX_adj, dataY_adj

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

def output_figure(gene_mean_x, gene_mean_y, adjusted_scores, monte_carlo_points, cer, start_x, start_y, grid_width, grid_height, units):
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

def csv_setup(fileName):
        # --- Load data ---

        # find number of genes
        numGene = -1
        numSample = -1
        with open(fileName, 'r') as file:
            csv_reader = csv.reader(file)
            for row in csv_reader:
                numGene+= 1
                if numSample == -1:
                    numSample = len(row) - 1
        file.close()

    
        # define profile variables
        if numGene < 1:
            print("error: could not read file")
        else:
            profile = np.zeros((numGene,numSample))
            sampleName = ['']*numSample
            geneName = ['']*numGene

            with open(fileName, 'r') as file:
                csv_reader = csv.reader(file)
                rowCount = 0
                for row in csv_reader:
                    rowCount+= 1
                    if rowCount == 1:
                        for i in range(numSample):
                            sampleName[i] = row[i+1]
                    else:
                        geneName[rowCount-2] = row[0]
                        profile[rowCount-2] = row[1:]

        controlInd = []
        treatmentInd = []

        for i in range(numSample):
            if sampleName[i] == 'control':
                controlInd.append(i)
            if sampleName[i] == 'treatment':
                treatmentInd.append(i)

        sampleName = np.array(sampleName)
        profile = np.array(profile)

        return profile[:,controlInd], profile[:,treatmentInd], geneName, sampleName



# Example usage
if __name__ == "__main__":
    # Test data
    data_x = np.random.rand(100,5)
    data_y = np.random.rand(100,5)
    mage(data_x, data_y)
