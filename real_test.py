# Matthew Beltran
# 1/20/2025

import MAGE
import numpy as np

def test():
    
    # --- Load profile from csv ---
    #control, treatment, geneName, sampleName = MAGE.csv_setup('MAGE Paper Examples\workspaces\gt3_NCBI.csv')
    control, treatment, geneName, sampleName = MAGE.csv_setup('MAGE Paper Examples\workspaces\mtor_NCBI.csv')
    
    # --- Run MAGE function ---
    #OS = MAGE.mage(control, treatment, output_plots=True, output_diags=False, saveFigs=False)

    # --- Calculate false discovery rate ---
    #FDR = MAGE.FDR(control, treatment, OS, output_plots=False, output_diags=False, saveFigs=True, contour_loops_max=5)
    #FDR = np.zeros(len(geneName))
    # --- Depth Test ---
    temp = MAGE.analyze_depth(control, treatment, (0.1, 0.3, 0.5, 0.7, 0.9))

    # --- Output to csv ---
    #np.savetxt('gt3_python_results.csv', np.c_[np.array(geneName),OS,FDR], fmt='%s', delimiter=',',header="Gene,OS,FDR")
    #np.savetxt('mtor_python_results.csv', np.c_[np.array(geneName),OS,FDR], fmt='%s', delimiter=',',header="Gene,OS,FDR")

if __name__ == "__main__":
    test()