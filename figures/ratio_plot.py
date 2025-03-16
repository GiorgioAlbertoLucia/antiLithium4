from ROOT import TRatioPlot, TCanvas, gStyle, TLegend
from torchic.core.histogram import load_hist, HistLoadInfo

if __name__ == '__main__':

    gStyle.SetOptStat(0)

    #h_CF23 = load_hist(HistLoadInfo('/home/galucia/antiLithium4/analysis/output/LHC23PbPb/studies.root', 'CorrelationAnti/hCorrelation_kstar'))
    #h_CF23.SetName('h_CF23')
    #h_CF23.SetMarkerColor(601)
    #h_CF23.SetMarkerStyle(20)
    #
    #h_CF24 = load_hist(HistLoadInfo('/home/galucia/antiLithium4/analysis/output/LHC24PbPb/studies.root', 'CorrelationAnti/hCorrelation_kstar'))
    #h_CF24.SetName('h_CF24')
    #h_CF24.SetTitle(';#it{k}* (GeV/#it{c});#it{C}(#it{k}*)')
    #h_CF24.SetMarkerColor(797)
    #h_CF24.SetMarkerStyle(20)

    h_CF = load_hist(HistLoadInfo('/home/galucia/antiLithium4/analysis/output/PbPb/studies.root', 'CorrelationAnti/hCorrelation_kstar'))
    h_CF.SetName('h_CF')
    #h_CF.SetMarkerColor(418)
    h_CF.SetMarkerColor(601)
    h_CF.SetMarkerStyle(20)

    h_CF_new = load_hist(HistLoadInfo('/home/galucia/antiLithium4/analysis/output/PbPb/studies_new.root', 'CorrelationAnti/hCorrelation_kstar'))
    h_CF_new.SetName('h_CF')
    h_CF_new.SetMarkerColor(797)
    h_CF_new.SetMarkerStyle(24)

    canvas = TCanvas('canvas', 'canvas', 800, 800)
    frame = canvas.DrawFrame(0.0, 0.0, 0.805, 1.5)
    #rp = TRatioPlot(h_CF24, h_CF23, 'divsym')
    rp = TRatioPlot(h_CF, h_CF_new, 'divsym')
    rp.SetH1DrawOpt('p e0 same')
    rp.SetH2DrawOpt('p e0 same')
    #rp.GetLowerPad().SetTitle('; #it{k}* (GeV/#it{c});Ratio (23/24)')
    rp.Draw('nogrid')

    rp.GetUpperPad().cd()
    legend = TLegend(0.5, 0.3, 0.8, 0.5)
    legend.SetBorderSize(0)
    #legend.AddEntry(h_CF23, '2023', 'p')
    #legend.AddEntry(h_CF24, '2024', 'p')
    legend.AddEntry(h_CF, '2023 + 2024', 'p')
    legend.AddEntry(h_CF_new, '(2023 + 2024) NEW', 'p')
    legend.Draw()
    #h_CF.Draw('p e0 same')

    #canvas.SaveAs('/home/galucia/antiLithium4/figures/24-02-2025/ratio_plot.pdf')
    canvas.SaveAs('/home/galucia/antiLithium4/figures/24-02-2025/ratio_plot_after_purity.pdf')

    
