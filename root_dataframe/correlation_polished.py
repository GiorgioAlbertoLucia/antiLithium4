
from ROOT import TFile, TH1F

infile_same = '/home/galucia/antiLithium4/root_dataframe/output/same_event_old.root'
infile_mixed = '/home/galucia/antiLithium4/root_dataframe/output/mixed_event_old.root'
outfile = TFile.Open('/home/galucia/antiLithium4/root_dataframe/output/correlation_old.root', 'RECREATE')

file_same = TFile.Open(infile_same)
file_mixed = TFile.Open(infile_mixed)



NORM_LOW_KSTAR = 0.2 # 0.25
NORM_HIGH_KSTAR = 0.4 # 0.75
NBINS_KSTAR = 40 # 20

for mode in ['', 'Matter', 'Antimatter']:

    for interval in ['']:

        outdir = outfile.mkdir(f'Correlation{mode}{interval}')

        # centrality = '0-10'
        h_same010 = file_same.Get(f'kstar{mode}/hKstar010{mode}{interval}')
        h_same010.SetDirectory(0)  # Detach from file to avoid issues with deletion
        h_same010.SetName('hSameEvent010')
        #h_same010.Rebin()


        h_mixed010 = file_mixed.Get(f'kstar{mode}/hKstar010{mode}{interval}')
        h_mixed010.SetDirectory(0)  # Detach from file to avoid issues with deletion
        h_mixed010.SetName('hMixedEvent010')
        #h_mixed010.Rebin()


        low_bin = h_same010.FindBin(NORM_LOW_KSTAR)
        high_bin = h_same010.FindBin(NORM_HIGH_KSTAR)
        normalization_factor = h_same010.Integral(low_bin, high_bin) / h_mixed010.Integral(low_bin, high_bin)
        h_mixed010.Scale(normalization_factor)

        h_corr010 = h_same010.Clone('hCorrelation010')
        h_corr010.Divide(h_mixed010)

        outdir.cd()
        h_same010.Write()
        h_mixed010.Write()
        h_corr010.Write('hCorrelationFull010')
        h_corr010_limited = TH1F('hCorrelation010Limited', 'hCorrelation010Limited', NBINS_KSTAR, 0, 0.4)
        for ibin in range(1, h_corr010_limited.GetNbinsX()+1):
            h_corr010_limited.SetBinContent(ibin, h_corr010.GetBinContent(ibin))
            h_corr010_limited.SetBinError(ibin, h_corr010.GetBinError(ibin))
        h_corr010_limited.Write('hCorrelation010')

        # centrality = '10-30'
        h_same1030 = file_same.Get(f'kstar{mode}/hKstar1030{mode}{interval}')
        h_same1030.SetDirectory(0)  # Detach from file to avoid issues with deletion
        h_same1030.SetName('hSameEvent1030')
        #h_same1030.Rebin()

        h_mixed1030 = file_mixed.Get(f'kstar{mode}/hKstar1030{mode}{interval}')
        h_mixed1030.SetDirectory(0)  # Detach from file to avoid issues with deletion
        h_mixed1030.SetName('hMixedEvent1030')
        #h_mixed1030.Rebin()

        low_bin = h_same1030.FindBin(NORM_LOW_KSTAR)
        high_bin = h_same1030.FindBin(NORM_HIGH_KSTAR)
        normalization_factor = h_same1030.Integral(low_bin, high_bin) / h_mixed1030.Integral(low_bin, high_bin)
        h_mixed1030.Scale(normalization_factor)

        h_corr1030 = h_same1030.Clone('hCorrelation1030')
        h_corr1030.Divide(h_mixed1030)

        outdir.cd()
        h_same1030.Write()
        h_mixed1030.Write()
        h_corr1030.Write('hCorrelationFull1030')
        h_corr1030_limited = TH1F('hCorrelation1030Limited', 'hCorrelation1030Limited', NBINS_KSTAR, 0, 0.4)
        for ibin in range(1, h_corr1030_limited.GetNbinsX()+1):
            h_corr1030_limited.SetBinContent(ibin, h_corr1030.GetBinContent(ibin))
            h_corr1030_limited.SetBinError(ibin, h_corr1030.GetBinError(ibin))
        h_corr1030_limited.Write('hCorrelation1030')


        # centrality = '30-50'
        h_same3050 = file_same.Get(f'kstar{mode}/hKstar3050{mode}{interval}')
        h_same3050.SetDirectory(0)  # Detach from file to avoid issues with deletion
        h_same3050.SetName('hSameEvent3050')
        #h_same3050.Rebin()

        h_mixed3050 = file_mixed.Get(f'kstar{mode}/hKstar3050{mode}{interval}')
        h_mixed3050.SetDirectory(0)  # Detach from file to avoid issues with deletion
        h_mixed3050.SetName('hMixedEvent3050')
        #h_mixed3050.Rebin()

        low_bin = h_same3050.FindBin(NORM_LOW_KSTAR)
        high_bin = h_same3050.FindBin(NORM_HIGH_KSTAR)
        normalization_factor = h_same3050.Integral(low_bin, high_bin) / h_mixed3050.Integral(low_bin, high_bin)
        h_mixed3050.Scale(normalization_factor)

        h_corr3050 = h_same3050.Clone('hCorrelation3050')
        h_corr3050.Divide(h_mixed3050)

        outdir.cd()
        h_same3050.Write()
        h_mixed3050.Write()
        h_corr3050.Write('hCorrelationFull3050')
        h_corr3050_limited = TH1F('hCorrelation3050Limited', 'hCorrelation3050Limited', NBINS_KSTAR, 0, 0.4)
        for ibin in range(1, h_corr3050_limited.GetNbinsX()+1):
            h_corr3050_limited.SetBinContent(ibin, h_corr3050.GetBinContent(ibin))
            h_corr3050_limited.SetBinError(ibin, h_corr3050.GetBinError(ibin))
        h_corr3050_limited.Write('hCorrelation3050')

        # centrality = '0-50%'
        h_same050 = h_same010.Clone('hSameEvent050')
        h_same050.Add(h_same1030)
        h_same050.Add(h_same3050)

        h_mixed050 = h_mixed010.Clone('hMixedEvent050')
        h_mixed050.Add(h_mixed1030)
        h_mixed050.Add(h_mixed3050)

        h_corr050 = h_same050.Clone('hCorrelation050')
        h_corr050.Divide(h_mixed050)

        outdir.cd()
        h_same050.Write()
        h_mixed050.Write()
        h_corr050.Write('hCorrelationFull050')
        h_corr050_limited = TH1F('hCorrelation050Limited', 'hCorrelation050Limited', NBINS_KSTAR, 0, 0.4)
        for ibin in range(1, h_corr050_limited.GetNbinsX()+1):
            h_corr050_limited.SetBinContent(ibin, h_corr050.GetBinContent(ibin))
            h_corr050_limited.SetBinError(ibin, h_corr050.GetBinError(ibin))
        h_corr050_limited.Write('hCorrelation050')

outfile.Close()