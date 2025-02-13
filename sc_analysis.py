# Matthew Beltran
# 2/5/2025

import MAGE
import numpy as np
import pickle

def test():
    
    # --- Load profile from pkl ---
    try:
        with open("sc_brain.pkl", "rb") as f:
            data = pickle.load(f)

        profile = data[0]
        geneName = data[1]
        sampleName = data[2]

        ECs = [i for i, s in enumerate(sampleName) if "EC" in s]
        Ms = [i for i, s in enumerate(sampleName) if "Mural" in s]

        # --- Run MAGE function ---
        #OS = MAGE.mage(profile[:,np.min(ECs):np.max(ECs)], profile[:,np.min(Ms):np.max(Ms)],
        #                output_plots=True, output_diags=False, saveFigs=False)
        #FDR = MAGE.FDR(profile[:,np.min(ECs):np.max(ECs)], profile[:,np.min(Ms):np.max(Ms)]
        #               , OS, output_plots=False, output_diags=False, saveFigs=True, contour_loops_max=5)
        temp = MAGE.analyze_samples(profile[:,np.min(ECs):np.max(ECs)], profile[:,np.min(Ms):np.max(Ms)],
                                   (0.001, 0.005, 0.01, 0.03, 0.05, 0.1, 0.3, 0.5), trials=5, saveData=True)
        #np.savetxt('mtor_python_results.csv', np.c_[np.array(geneName),OS,FDR], fmt='%s', delimiter=',',header="Gene,OS,FDR")


    except EOFError:
        print('Data not found...')

    # --- Output to csv ---
    #np.savetxt('gt3_python_results.csv', np.c_[np.array(geneName),OS,FDR], fmt='%s', delimiter=',',header="Gene,OS,FDR")
    #np.savetxt('mtor_python_results.csv', np.c_[np.array(geneName),OS,FDR], fmt='%s', delimiter=',',header="Gene,OS,FDR")

if __name__ == "__main__":
    test()