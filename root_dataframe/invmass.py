from ROOT import TFile, TH1F
import numpy as np


def pull(h1, h2):
    """
    Calculate the pull between two histograms.
    """
    if h1.GetNbinsX() != h2.GetNbinsX():
        print(f'{h1.GetNbinsX()=}, {h2.GetNbinsX()=}')
        raise ValueError("Histograms must have the same number of bins")

    pull_hist = h1.Clone("hPull")
    pull_hist.Reset()

    for bin in range(1, h1.GetNbinsX() + 1):
        content1 = h1.GetBinContent(bin)
        content2 = h2.GetBinContent(bin)
        error1 = h1.GetBinError(bin)
        error2 = h2.GetBinError(bin)

        if error1 > 0 and error2 > 0:
            pull_value = (content1 - content2) / np.sqrt(error1**2 + error2**2)
            pull_hist.SetBinContent(bin, pull_value)
            pull_hist.SetBinError(bin, 1)  # Pulls are typically not assigned errors

    return pull_hist

def subtraction(h1, h2):
    """
    Subtract histogram h2 from h1.
    """
    if h1.GetNbinsX() != h2.GetNbinsX():
        print(f'{h1.GetNbinsX()=}, {h2.GetNbinsX()=}')
        raise ValueError("Histograms must have the same number of bins")

    sub_hist = h1.Clone("hSubtraction")
    sub_hist.Reset()

    for bin in range(1, h1.GetNbinsX() + 1):
        content1 = h1.GetBinContent(bin)
        content2 = h2.GetBinContent(bin)
        error1 = h1.GetBinError(bin)
        error2 = h2.GetBinError(bin)
        sub_hist.SetBinContent(bin, content1 - content2)
        sub_hist.SetBinError(bin, np.sqrt(error1**2 + error2**2))

    return sub_hist

infile_same = '/home/galucia/antiLithium4/root_dataframe/output/same_event.root'
infile_mixed = '/home/galucia/antiLithium4/root_dataframe/output/mixed_event.root'
outfile = TFile.Open('/home/galucia/antiLithium4/root_dataframe/output/invariant_mass.root', 'RECREATE')

file_same = TFile.Open(infile_same)
file_mixed = TFile.Open(infile_mixed)

for mode in ['', 'Matter', 'Antimatter']:

    outdir = outfile.mkdir(f'InvariantMass{mode}')

    h_same = file_same.Get(f'InvariantMass{mode}/hInvariantMass{mode}')
    h_same.SetDirectory(0)  # Detach from file to avoid issues with deletion
    h_same.SetName('hSameEvent')
    h_same.Rebin()

    h_mixed = file_mixed.Get(f'InvariantMass{mode}/hInvariantMass{mode}')
    h_mixed.SetDirectory(0)  # Detach from file to avoid issues with deletion
    h_mixed.SetName('hMixedEvent')
    h_mixed.Rebin()

    low_bin = h_same.FindBin(3.78)
    high_bin = h_same.FindBin(3.85)
    normalization_factor = h_same.Integral(low_bin, high_bin) / h_mixed.Integral(low_bin, high_bin)
    h_mixed.Scale(normalization_factor)

    h_pull = pull(h_same, h_mixed)
    h_subtraction = subtraction(h_same, h_mixed)

    outdir.cd()
    h_same.Write()
    h_mixed.Write()
    h_pull.Write()
    h_subtraction.Write()