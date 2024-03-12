import numpy as np
import pandas as pd
import uproot
from ROOT import TFile, TH1F

HIST = {
    #'fSignedPtHe3': ['fSignedPtHe3', 100, -5, 5],
    #'fSignedPtPr': ['fSignedPtPr', 100, -5, 5],
    'fIsMatter': ['fIsMatter', 2, 0, 2],
    'fHe3NSigmaTPC': ['fHe3NSigmaTPC', 100, -5, 5],
    'fPrNSigmaTPC': ['fPrNSigmaTPC', 100, -5, 5],
    'fHe3MassTOF': ['fHe3MassTOF', 100, 2.9, 3.2],
    'fPrMassTOF': ['fPrMassTOF', 100, 0.7, 1.1],
    'fMassMC': ['fMassMC', 100, 0, 5],
    'fSignedPtMC': ['fSignedPtMC', 100, -5, 5]
}

def visualizeDF(uprootFile, dfID, visualHISTs):

    tree = uprootFile[f'DF_{dfID}']['O2lithium4tablemc']
    #tree.arrays('fSignedPtHe3', aliases={'fSignedPtHe3': 'fPtHe3 * ((fIsMatter == 0) - (fIsMatter == 1))'})
    #tree.arrays('fSignedPtPr', aliases={'fSignedPtPr': 'fPtPr * ((fIsMatter == 0) - (fIsMatter == 1))'})
    
    #df = tree.arrays(['fSignedPtHe3', 'fSignedPtPr', 'fHe3NSigmaTPC', 'fPrNSigmaTPC', 'fHe3MassTOF', 'fPrMassTOF', 'fSignedPtMC'], library='pd')
    df = tree.arrays(['fIsMatter', 'fHe3NSigmaTPC', 'fPrNSigmaTPC', 'fHe3MassTOF', 'fPrMassTOF', 'fMassMC', 'fSignedPtMC'], library='pd')
    print(df.describe())

    for key in HIST:    visualHISTs[key].FillN(len(df[key]), np.asarray(df[key].values, dtype=float), np.ones(df[key].size, dtype=float))


    return visualHISTs


def main(inFilePath, visualPath):

    inFile = uproot.open(inFilePath)
    visualFile = TFile(visualPath, "RECREATE")
    visualHISTs = {key: TH1F(key, HIST[key][0], HIST[key][1], HIST[key][2], HIST[key][3]) for key in HIST} 

    visualHISTs = visualizeDF(inFile, 1, visualHISTs)
    visualHISTlist = []
    
    
    for key in visualHISTs:     visualHISTlist.append(visualHISTs[key])
    visualFile.cd()
    for hist in visualHISTlist: 
        hist.Write()
        hist.Print('base')
        
    visualFile.Close()



if __name__ == "__main__":

    #main("/home/galucia/antiLithium4/MCWorkflowLauncher/AO2D_lit_mc.root", "/home/galucia/antiLithium4/MCanalysis/output/visual.root")
    main("/Users/glucia/Projects/ALICE/antiLithium4/MCWorkflowLauncher/AO2D_lit_mc.root", 
         "/Users/glucia/Projects/ALICE/antiLithium4/MCanalysis/output/visual.root")
    
