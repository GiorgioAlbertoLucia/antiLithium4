'''
    Script to convert the histograms in a file from MeV to GeV.
'''

from ROOT import TH1F, TFile

def convert_MeV_to_GeV(hist: TH1F):
    '''
        Convert the histogram from MeV to GeV.
    '''

    hist_name = hist.GetName()
    hist.SetName(f'{hist_name}_to_convert')
    converted_hist = TH1F(hist_name, hist.GetTitle(), hist.GetNbinsX(), hist.GetXaxis().GetXmin() / 1000., hist.GetXaxis().GetXmax() / 1000.)

    for ibin in range(1, hist.GetNbinsX() + 1):
        converted_hist.SetBinContent(ibin, hist.GetBinContent(ibin))
        converted_hist.SetBinError(ibin, hist.GetBinError(ibin))
    
    return converted_hist

if __name__ == '__main__':

    infiles = ['CATS_CF_LS_5p00fm.root',
               'CATS_CF_LS_4p45fm.root',
               'CATS_CF_LS_4p16fm.root',
               'CATS_CF_LS_3p19fm.root'] 
    hist_names = ['hHe3_p_Coul_CF', 'hHe3_p_Coul_InvMass']

    for infile in infiles:

        read_file = TFile.Open(infile, 'READ')
        write_file = TFile.Open(infile.replace('.root', '_converted.root'), 'RECREATE')

        for hist_name in hist_names:
            hist = read_file.Get(hist_name)
            converted_hist = convert_MeV_to_GeV(hist)
            write_file.cd()
            converted_hist.Write()
        
        write_file.Close()
        read_file.Close()