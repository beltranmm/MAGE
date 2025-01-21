# Matthew Beltran
# 1/20/2024

import MAGE
import numpy as np
import csv

def test():
    
    # --- Load data ---

    # find number of genes
    numGene = -1
    numSample = -1
    with open('MAGE Paper Examples\workspaces\mtor_NCBI.csv', 'r') as file:
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

        with open('MAGE Paper Examples\workspaces\mtor_NCBI.csv', 'r') as file:
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
    
    # --- Run MAGE function ---
    print("Running gt3 test...")
    OS = MAGE.mage(profile[:,controlInd], profile[:,treatmentInd])
    #FDR = MAGE.FDR(data_x, data_y, OS)


if __name__ == "__main__":
    test()