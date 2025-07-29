from ROOT import TFile
from torchic import HistLoadInfo
from torchic.core.histogram import load_hist

if __name__ == '__main__':

    infile_paths = ['/home/galucia/antiLithium4/root_dataframe/output/mixed_event_old_23.root',
                    '/home/galucia/antiLithium4/root_dataframe/output/mixed_event_old_24.root']
    outfile = TFile.Open('/home/galucia/antiLithium4/root_dataframe/output/mixed_event_old.root', 'RECREATE')
    plots = [
                'kstar/hKstar010', 'kstar/hKstar1030',
                'kstar/hKstar3050',
                'kstarMatter/hKstar010Matter', 'kstarMatter/hKstar1030Matter',
                'kstarMatter/hKstar3050Matter',
                'kstarAntimatter/hKstar010Antimatter', 'kstarAntimatter/hKstar1030Antimatter',
                'kstarAntimatter/hKstar3050Antimatter',
            ]

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
        