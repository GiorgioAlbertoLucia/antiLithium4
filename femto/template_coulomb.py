'''
    Produce a template fit for the Coulomb interaction for a p-He3 system
'''

from torchic.core.histogram import load_hist, HistLoadInfo
from ROOT import RooRealVar, RooDataHist, RooHistPdf
from ROOT import TFile


if __name__ == '__main__':

    infile_path = '/home/galucia/antiLithium4/analysis/output/CATS/CATS_cent0_10_new.root'
    hist = load_hist(HistLoadInfo(infile_path, 'hHe3_p_Coul_CF_LS'))

    kstar = RooRealVar('kstar', 'kstar', 0, 0.8)
    data_hist = RooDataHist('data_hist', 'data_hist', [kstar], Import=hist)
    pdf = RooHistPdf('pdf', 'pdf', {kstar}, data_hist, 0)

    frame = kstar.frame()
    data_hist.plotOn(frame)
    pdf.plotOn(frame)

    outfile = TFile.Open('/home/galucia/antiLithium4/femto/output/template_coulomb.root', 'recreate')
    frame.Write()
    outfile.Close()
