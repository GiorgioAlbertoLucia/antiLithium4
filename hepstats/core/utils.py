import numpy as np
from ROOT import TH1F
from hist import Hist
import hist
import zfit
import boost_histogram as bh
from copy import deepcopy

def hist_to_tensor(hist: TH1F, xmin, xmax):
    '''
        Converts a TH1 histogram to a tensor.
    '''

    bin_min = hist.FindBin(xmin)
    bin_max = hist.FindBin(xmax)
    nbins = bin_max - bin_min + 1
    xedges = np.zeros(nbins+1)
    ys, yerrs = np.zeros(nbins), np.zeros(nbins)

    for ibin in range(nbins):
        xedges[ibin] = hist.GetXaxis().GetBinLowEdge(bin_min + ibin)
        ys[ibin] = hist.GetBinContent(bin_min + ibin)
        yerrs[ibin] = hist.GetBinError(bin_min + ibin)
    xedges[ibin+1] = hist.GetXaxis().GetBinUpEdge(bin_max)

    return xedges, ys, yerrs

def pdf_to_hist(pdf, xmin:float, xmax:float, title:str=''):
    '''
        Converts a zfit pdf to a histogram.
    '''
    
    nbins = 1000
    xs = np.linspace(xmin, xmax, nbins).reshape(-1, 1)
    ys = pdf.pdf(xs).numpy().flatten()
    
    hist = Hist.new.Reg(nbins, xmin, xmax, name=f'{pdf.name}_hist').Double()
    hist.fill(xs.flatten(), weight=ys)

    return hist

def pdf_to_numpy(pdf, obs:zfit.Space, nbins:int=1000):
    '''
        Converts a zfit pdf to a histogram.
    '''
    
    xs = np.linspace(obs.lower, obs.upper, nbins).reshape(-1, 1)
    ys = pdf.pdf(xs, norm_range=obs).numpy().flatten()
    
    return xs, ys
    
def space_from_hist(hist: Hist, name:str=''):
    '''
        Create a space from a histogram binning.
    '''

    bin_edges = [hist.axes.edges[0][0], hist.axes.edges[0][-1]]
    nbins = len(hist.axes.edges[0]) - 1
    print(f'{bin_edges=}')
    binning = zfit.binned.RegularBinning(nbins, bin_edges[0], bin_edges[1], name=name)
    return zfit.Space(name, binning=binning)

def ratio_hist(hist1: Hist, hist2: Hist):
    '''
        Calculate the ratio of two histograms.
    '''

    h_ratio = Hist(
        *hist1.axes,
        storage=bh.storage.Weight()
    )

    vals = hist1.values() / hist2.values()
    var_num = hist1.variances()
    var_den = hist2.variances()
    variances = (vals**2) * \
          ((var_num / hist1.values()**2) + (var_den / hist2.values()**2))
    
    view = h_ratio.view()
    view.value[:]    = vals
    view.variance[:] = variances

    return h_ratio

def match_hist_binning(source_hist:Hist, target_hist:Hist):
    '''
        Match the binning of hist1 to hist2.
    '''

    target_edges = target_hist.axes[0].edges
    rebinned_hist = Hist(hist.axis.Variable(target_edges, name=source_hist.axes[0].name))

    source_edges = source_hist.axes[0].edges
    source_centers = (source_edges[:-1] + source_edges[1:]) / 2
    source_values = source_hist.values()

    data = np.repeat(source_centers, source_values.astype(int))
    rebinned_hist.fill(data)

    return rebinned_hist