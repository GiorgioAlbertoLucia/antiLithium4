from ROOT import TFile
from torchic import HistLoadInfo
from torchic.core.histogram import load_hist

if __name__ == '__main__':

    infile_paths = ['/home/galucia/antiLithium4/analysis/output/LHC23PbPb/event_mixing_visual_noH3_selectionsPr.root',
                    '/home/galucia/antiLithium4/analysis/output/LHC24PbPb/event_mixing_visual_noH3_selectionsPr.root']
    outfile = TFile.Open('/home/galucia/antiLithium4/analysis/output/PbPb/event_mixing_visual_noH3_selectionsPr.root', 'RECREATE')
    plots = ['Correlations/fKstarCentralityAnti', 'Correlations/fKstarCentralityMatter', 
             'InvMass/InvMassCentralityMatter', 'InvMass/InvMassCentralityAnti']

    for plot in plots:
        dir_name = plot.split('/')[0]
        dir = outfile.GetDirectory(dir_name)
        if not dir:
            dir = outfile.mkdir(dir_name)
        
        hist = None
        for infile_path in infile_paths:
            if not hist:
                hist = load_hist(HistLoadInfo(infile_path, plot))
            else:
                tmp_hist = load_hist(HistLoadInfo(infile_path, plot))
                hist.Add(tmp_hist)
        dir.cd()
        hist.Write()
        