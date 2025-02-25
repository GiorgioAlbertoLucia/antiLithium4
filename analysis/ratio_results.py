'''
    Ratio of the results of LHC24PbPb and LHC23PbPb
'''

from torchic import HistLoadInfo
from torchic.core.histogram import load_hist
from ROOT import TFile

def ratio_hist(h1, h2):
    '''Ratio of two histograms'''
    hist = h1.Clone(f'ratio_{h1.GetName()}')
    hist.Divide(h2)
    return hist

if __name__ == '__main__':

    infile23_path = '/home/galucia/antiLithium4/analysis/output/LHC23PbPb/studies_noH3.root'
    infile24_path = '/home/galucia/antiLithium4/analysis/output/LHC24PbPb/studies_noH3.root'
    outfile = TFile.Open('/home/galucia/antiLithium4/analysis/output/PbPb/ratio_studies_noH3.root', 'RECREATE')

    plots = ['CorrelationMatter/hCorrelation_kstar', 
             'CorrelationMatter/hCorrelation_kstar_cent0.0_10.0',
             'CorrelationMatter/hCorrelation_kstar_cent10.0_30.0',
             'CorrelationMatter/hCorrelation_kstar_cent30.0_50.0',
             'CorrelationAnti/hCorrelation_kstar', 
             'CorrelationAnti/hCorrelation_kstar_cent0.0_10.0',
             'CorrelationAnti/hCorrelation_kstar_cent10.0_30.0',
             'CorrelationAnti/hCorrelation_kstar_cent30.0_50.0',
            ]

    for plot in plots:
        h23 = load_hist(HistLoadInfo(infile23_path, plot))
        h23.SetName(f'{h23.GetName()}_{plot.split("/")[0][11:]}')
        h24 = load_hist(HistLoadInfo(infile24_path, plot))
        ratio = ratio_hist(h23, h24)
        outfile.cd()
        ratio.Write()

    outfile.Close()