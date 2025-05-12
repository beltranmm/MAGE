# Matthew Beltran
# 1/20/2025

import MAGE
import numpy as np

import matplotlib.pyplot as plt

def test():
    
    # --- Load profile from csv ---
    #control, treatment, geneName, sampleName = MAGE.csv_setup('MAGE Paper Examples\workspaces\gt3_NCBI.csv')
    control, treatment, geneName, sampleName = MAGE.csv_setup('MAGE Paper Examples\workspaces\mtor_NCBI.csv')
    
    #Pabpn1 = np.where(np.asarray(geneName) == 'Pabpn1')
    #print(geneName[Pabpn1[0][0]])
    #plt.figure(figsize=(10, 6))
    #plt.scatter(control[Pabpn1[0][0],:], treatment[Pabpn1[0][0],:],
    #             c='g', s=5, edgecolor='g', alpha=0.5, label='MC Sampling')
    #plt.show()

    # --- Run MAGE function ---
    #OS = MAGE.mage(control, treatment, remove_high_low_expr=False, dist_adj=False, output_plots=True, output_diags=False, saveFigs=True)
    #OS = MAGE.mage(control, treatment, contour_loops_max=10, num_starting_contours=200, nonnegative=True, remove_high_low_expr=True, dist_adj=True, output_plots=True, output_diags=False, saveFigs=False)

    # --- Calculate false discovery rate ---
    #FDR = MAGE.FDR(control, treatment, OS, remove_high_low_expr=False, dist_adj=False, output_plots=True, output_diags=False, saveFigs=False)
    #FDR = np.zeros(len(geneName))
    # --- Depth Test ---
    temp = MAGE.analyze_depth(control, treatment, (0.5,1), save_data=True)
    #temp = MAGE.analyze_samples(control, treatment, (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9))

    # --- Output to csv ---
    #np.savetxt('gt3_python_results.csv', np.c_[np.array(geneName),OS,FDR], fmt='%s', delimiter=',',header="Gene,OS,FDR")
    #np.savetxt('mtor_python_results.csv', np.c_[np.array(geneName),OS,FDR], fmt='%s', delimiter=',',header="Gene,OS,FDR")

if __name__ == "__main__":
    test()