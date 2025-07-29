'''
    Simple script to perform corrected subtraction of the invariant mass histograms
'''

import numpy as np
import uproot
from ROOT import TH1F, TFile
from torchic.core.histogram import HistLoadInfo, load_hist

def correct_mixed_event(h_mixed_event:TH1F, h_correction:TH1F):

    for ibin in range(1, h_mixed_event.GetNbinsX()+1):

        xvalue = h_mixed_event.GetBinCenter(ibin)
        yvalue = h_mixed_event.GetBinContent(ibin)

        jbin = h_correction.FindBin(xvalue)
        ycorrection = h_correction.GetBinContent(jbin)

        h_mixed_event.SetBinContent(ibin, yvalue*ycorrection)
        h_mixed_event.SetBinError(ibin, np.sqrt(yvalue*ycorrection))

    return h_mixed_event

def normalise_mixed_event(h_same_event:TH1F, h_mixed_event: TH1F, xlow: float, xup: float):

    binlow_same = h_same_event.FindBin(xlow)
    binup_same = h_same_event.FindBin(xup)
    integral_same_event = h_same_event.Integral(binlow_same, binup_same)

    binlow_mixed = h_mixed_event.FindBin(xlow)
    binup_mixed = h_mixed_event.FindBin(xup)
    integral_mixed_event = h_mixed_event.Integral(binlow_mixed, binup_mixed)

    normalisation = integral_same_event / integral_mixed_event

    h_mixed_event.Scale(normalisation)
    return h_mixed_event

def mev_to_gev(hist:TH1F):
    '''
        Convert the x-axis of the histogram from MeV to GeV.
    '''
    
    xmin_gev = hist.GetXaxis().GetXmin() / 1000.    
    xmax_gev = hist.GetXaxis().GetXmax() / 1000.
    nbins = hist.GetNbinsX()

    hist_gev = TH1F(hist.GetName()+'_gev', hist.GetTitle(), nbins, xmin_gev, xmax_gev)
    for ibin in range(1, nbins+1):
        hist_gev.SetBinContent(ibin, hist.GetBinContent(ibin))

    return hist_gev

def ratio_hist(h1: TH1F, h2: TH1F):
    '''
        Calculate the ratio of two histograms.
    '''
    
    if h1.GetNbinsX() != h2.GetNbinsX():
        raise ValueError("Histograms must have the same number of bins")

    h_ratio = h1.Clone(h1.GetName()+'_ratio')
    for ibin in range(1, h1.GetNbinsX()+1):
        if h2.GetBinContent(ibin) != 0:
            ratio = h1.GetBinContent(ibin) / h2.GetBinContent(ibin)
            h_ratio.SetBinContent(ibin, ratio)
            h_ratio.SetBinError(ibin, ratio*np.sqrt(1/h1.GetBinContent(ibin) + 1/h2.GetBinContent(ibin)))
        else:
            h_ratio.SetBinContent(ibin, 0)
            h_ratio.SetBinError(ibin, 0)

    return h_ratio


if __name__ == '__main__':

    h_same_event = load_hist(HistLoadInfo(
        '/home/galucia/antiLithium4/analysis/output/LHC24PbPb/data_visual_selectionsPr.root',
        'InvMass/InvMassAntiLi'))
    h_same_event.Rebin(2)

    h_mixed_event = load_hist(HistLoadInfo(
        '/home/galucia/antiLithium4/analysis/output/LHC24PbPb/event_mixing_visual_selectionsPr.root',
        'InvMass/InvMassAntiLi'))
    h_mixed_event.Rebin(2)

    h_correction = load_hist(HistLoadInfo(
        '/home/galucia/antiLithium4/analysis/output/CATS/CATS_CF_LS.root',
        'hHe3_p_Coul_InvMass'
    ))
    h_correction = mev_to_gev(h_correction)

    outfile = TFile.Open('output/subtraction.root', 'recreate')

    h_mixed_event.Write('h_mixed_event_before_correction')
    h_mixed_event = correct_mixed_event(h_mixed_event, h_correction)
    h_mixed_event.Write('h_mixed_event_before_normalisation')
    h_mixed_event = normalise_mixed_event(h_same_event, h_mixed_event, 3.78, 3.84)
    
    h_subtraction = h_same_event.Clone('h_subtraction')
    h_subtraction.Add(h_mixed_event, -1.)
    h_ratio = ratio_hist(h_same_event, h_mixed_event)


    h_same_event.Write('h_same_event')
    h_mixed_event.Write('h_mixed_event')
    h_correction.Write('h_correction')
    h_subtraction.Write()
    h_ratio.Write('h_ratio')

    outfile.Close()
