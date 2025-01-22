'''
    Compare the outputs of the two event mixing methods
'''

import numpy as np
from ROOT import TFile, TCanvas


if __name__ == '__main__':

    em_file_grid = TFile.Open('/home/galucia/antiLithium4/analysis/output/LHC23PbPb/event_mixing_visual_grid_selectionsPr.root')
    em_file = TFile.Open('/home/galucia/antiLithium4/analysis/output/LHC23PbPb/event_mixing_visual_LSonly_selectionsPr.root')
    outfile = TFile.Open('/home/galucia/antiLithium4/analysis/output/LHC23PbPb/event_mixing_visual_comparison.root', 'RECREATE')

    dirs = ["QA", "InvMass", "Kinematics", "ITS", "TPC", "TOF", "DCA", "PID", "Centrality", "Efficiency", "Correlations", "SelfCorrelations"]
    for dir in dirs:
        print(f'\nProcessing {dir}')
        keys = em_file[dir].GetListOfKeys()
        for key in keys:
            print(f'Processing {key.GetName()}')

            hist_grid = em_file_grid[dir].Get(key.GetName())
            hist = em_file[dir].Get(key.GetName())
            hist_ratio = hist.Clone(f'{hist.GetName()}_ratio')
            
            if 'TH1' not in str(type(hist_grid)) or 'TH1' not in str(type(hist)):
                continue
            
            integral_grid = hist_grid.Integral(1, hist_grid.GetNbinsX()+1)
            if 'Kstar' in hist_grid.GetName():
                integral_grid = hist_grid.Integral(0, 160)
            if 'InvMass' in hist_grid.GetName():
                integral_grid = hist_grid.Integral(0, 162)
            integral = hist.Integral()
            if integral_grid > 0:
                hist_grid.Scale(integral/integral_grid)

            hist_ratio.Reset()
            for ibin in range(1, hist.GetNbinsX()+1):
                if hist_grid.GetBinContent(ibin) > 0:
                    hist_ratio.SetBinContent(ibin, hist.GetBinContent(ibin)/hist_grid.GetBinContent(ibin))
                    error = np.sqrt((hist.GetBinError(ibin)/hist_grid.GetBinContent(ibin))**2 + (hist_grid.GetBinError(ibin)*hist_ratio.GetBinContent(ibin)/hist_grid.GetBinContent(ibin))**2)
                    hist_ratio.SetBinError(ibin, error)
            
            outfile.cd()
            hist_ratio.Write()
    
