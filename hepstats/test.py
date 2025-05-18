'''
    Test the drawing of histograms and the fitting
'''

from torchic import Dataset, AxisSpec, HistLoadInfo
from torchic.core.histogram import load_hist

import uproot
import matplotlib.pyplot as plt
import mplhep
import zfit
from hist import Hist
from core.ehist import EHist

if __name__ == '__main__':

    #dataset = Dataset.from_root('/home/galucia/antiLithium4/analysis/output/LHC24PbPb/data.root', tree_name='outTree')
    #print(dataset.columns)
    #dataset.query('fIsMatter == 0 and fKstar < 0.8', inplace=True)
    
    #h_bkg = load_hist(HistLoadInfo(f'/home/galucia/antiLithium4/analysis/output/CATS/CATS_new.root', 
    #                                                       'hHe3_p_Coul_CF_LS'))
    file = uproot.open('/home/galucia/antiLithium4/analysis/output/CATS/CATS_new.root')
    h_bkg = file['hHe3_p_Coul_CF_LS'].to_hist()
    
    kstar = zfit.Space('kstar', limits=(0., 0.3))
    #hb_bkg = EHist.from_root(h_bkg, name='hb_bkg')
    #hb_bkg = EHist(Hist.new.Reg(h_bkg.GetNbinsX(), h_bkg.GetXaxis().GetXmin(), 
    #                       h_bkg.GetXaxis().GetXmax(), name='hb_bkg').Double())
    #hb_bkg.fill_from_root(h_bkg)
    
    
    #bkg_counts = zfit.Parameter('bkg_counts', 1., 0., 1e6)
    bkg_pdf = zfit.pdf.HistogramPDF(h_bkg, norm=None, name='bkg_pdf')

    # plot 
    mplhep.histplot(h_bkg, yerr=h_bkg.variances, label='h_bkg', histtype='step', color='blue')
    mplhep.histplot(bkg_pdf.to_hist(), label='bkg_pdf', histtype='step', color='red')
    plt.xlabel(r'\it{k}^{*} (GeV/\it{c})')
    plt.legend()
    plt.savefig(f'output/bkg_pdf.pdf')

    with uproot.recreate('output/test.root') as file:
        file['h_bkg'] = h_bkg
        file['bkg_pdf'] = bkg_pdf.to_hist()
