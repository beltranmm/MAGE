# Matthew Beltran
# 1/20/2024

import MAGE

def test():
    
    # --- Load profile from csv ---
    #control, treatment = MAGE.csv_setup('MAGE Paper Examples\workspaces\gt3_NCBI.csv')
    control, treatment = MAGE.csv_setup('MAGE Paper Examples\workspaces\mtor_NCBI.csv')
    
    # --- Run MAGE function ---
    OS = MAGE.mage(control, treatment, output_plots=False, output_diags=False, saveFigs=True)

    # --- Calculate false discovery rate ---
    FDR = MAGE.FDR(control, treatment, OS, output_plots=False, output_diags=False, saveFigs=True)


if __name__ == "__main__":
    test()