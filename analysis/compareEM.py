'''
    Compare the outputs of the two event mixing methods
'''

import numpy as np
from ROOT import TFile, TCanvas

def ratio_comparison():
    pass

def superimpose_comparison(outfile: TFile):
    
    se_file = TFile.Open('/Users/glucia/Projects/ALICE/antiLithium4/analysis/output/LHC24PbPb/data_visual_selectionsPr.root')
    me_file = TFile.Open('/Users/glucia/Projects/ALICE/antiLithium4/analysis/output/LHC24PbPb/event_mixing_visual_selectionsPr.root')
    #me_file_grid = TFile.Open('/Users/glucia/Projects/ALICE/antiLithium4/analysis/output/LHC24PbPb/event_mixing_visual_grid_selectionsPr.root')
    me_file_local_ar = TFile.Open('/Users/glucia/Projects/ALICE/data/lithium/mixing/LHC24ar_pass1_mixing_small.root')
    me_file_local_as = TFile.Open('/Users/glucia/Projects/ALICE/data/lithium/mixing/LHC24as_pass1_mixing_small.root')
    plots = ['Kinematics/PtAntiHe3', 'Kinematics/PtAntiHad', 'DCA/DCAxyHe3', 'DCA/DCAxyHad']
    
    for plot in plots:
        canvas = TCanvas(plot.split('/')[-1], plot.split('/')[-1])

        se_hist = se_file.Get(plot)
        se_hist.Scale(1/se_hist.Integral())
        se_hist.SetLineColor(2)
        se_hist.SetTitle('same-event;#it{p}_{T} (GeV/c);')
        
        me_hist = me_file.Get(plot)
        me_hist.Scale(1/me_hist.Integral())
        me_hist.SetLineColor(4) 
        me_hist.SetTitle('mixed-event-local;#it{p}_{T} (GeV/c);')

        if plot == 'Kinematics/PtAntiHad':
            se_int = se_hist.Integral(0, se_hist.FindBin(0))
            me_hist.Scale(se_int)

        if plot == 'Kinematics/PtAntiHe3':
            me_hist_local_ar = me_file_local_ar.Get('hHe3Unique')
            me_hist_local_as = me_file_local_as.Get('hHe3Unique')
            me_hist_local = me_hist_local_ar.Clone('hHe3Unique_local')
            me_hist_local.Add(me_hist_local_as)
            me_hist_local.Rebin()
            me_hist_local.Scale(1/me_hist_local.Integral())
            me_hist_local.SetLineColor(6)
            me_hist_local.SetTitle('mixed-event-local-unique;#it{p}_{T} (GeV/c);')

        canvas.cd()
        me_hist.Draw('hist same')
        se_hist.Draw('hist same')
        if plot == 'Kinematics/PtAntiHe3':
            me_hist_local.Draw('hist same')
        
        outfile.cd()
        canvas.BuildLegend()
        canvas.Write()




if __name__ == '__main__':

    #em_file_grid = TFile.Open('/home/galucia/antiLithium4/analysis/output/LHC23PbPb/event_mixing_visual_grid_selectionsPr.root')
    #em_file = TFile.Open('/home/galucia/antiLithium4/analysis/output/LHC23PbPb/event_mixing_visual_selectionsPr.root')
    #outfile = TFile.Open('/home/galucia/antiLithium4/analysis/output/LHC23PbPb/event_mixing_visual_comparison.root', 'RECREATE')
#
    #dirs = ["QA", "InvMass", "Kinematics", "ITS", "TPC", "TOF", "DCA", "PID", "Centrality", "Efficiency", "Correlations", "SelfCorrelations"]
    #for dir in dirs:
    #    print(f'\nProcessing {dir}')
    #    keys = em_file[dir].GetListOfKeys()
    #    for key in keys:
    #        print(f'Processing {key.GetName()}')
#
    #        hist_grid = em_file_grid[dir].Get(key.GetName())
    #        hist = em_file[dir].Get(key.GetName())
    #        hist_ratio = hist.Clone(f'{hist.GetName()}_ratio')
    #        
    #        if 'TH1' not in str(type(hist_grid)) or 'TH1' not in str(type(hist)):
    #            continue
    #        
    #        integral_grid = hist_grid.Integral()
    #        if 'Kstar' in hist_grid.GetName():
    #            integral_grid = hist_grid.Integral(0, 160)
    #        if 'InvMass' in hist_grid.GetName():
    #            integral_grid = hist_grid.Integral(0, 162)
    #        integral = hist.Integral()
    #        if integral_grid > 0:
    #            hist_grid.Scale(integral/integral_grid)
#
    #        hist_ratio.Reset()
    #        for ibin in range(1, hist.GetNbinsX()+1):
    #            if hist_grid.GetBinContent(ibin) > 0:
    #                hist_ratio.SetBinContent(ibin, hist.GetBinContent(ibin)/hist_grid.GetBinContent(ibin))
    #                error = np.sqrt((hist.GetBinError(ibin)/hist_grid.GetBinContent(ibin))**2 + (hist_grid.GetBinError(ibin)*hist_ratio.GetBinContent(ibin)/hist_grid.GetBinContent(ibin))**2)
    #                hist_ratio.SetBinError(ibin, error)
    #        
    #        outfile.cd()
    #        hist_ratio.Write()

    outfile = TFile.Open('/Users/glucia/Projects/ALICE/antiLithium4/analysis/output/LHC24PbPb/compareEM.root', 'RECREATE')
    superimpose_comparison(outfile)
    outfile.Close()
