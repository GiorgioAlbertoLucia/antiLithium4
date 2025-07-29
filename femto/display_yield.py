from ROOT import TFile, TCanvas, TPaveText, TLegend, gStyle, kRed, kBlue, kBlack, TVectorD
from torchic.core.histogram import load_hist, HistLoadInfo

def draw_same_event_li():

    h_same_event_norm = load_hist(
        HistLoadInfo('/home/galucia/antiLithium4/femto/output/AntiLithium4FitCFCATS.root',
        'dir/h_same_event_norm1.0')
    )

    gStyle.SetOptStat(0)
    canvas = TCanvas('canvas', 'canvas', 800, 600)
    #canvas.SetLogy()

    h_same_event_norm.SetMarkerColor(797)
    h_same_event_norm.SetMarkerStyle(20)
    h_same_event_norm.SetLineColor(797)
    h_same_event_norm.SetLineWidth(2)
    h_same_event_norm.Draw('hist')
    h_same_event_norm.SetTitle(';#it{k}* (GeV/#it{c}); Counts')

    text = TPaveText(0.55, 0.5, 0.85, 0.75, 'NDC')
    text.SetFillColor(0)
    text.SetBorderSize(0)
    text.AddText('This work')
    text.AddText('#bf{ALICE, Run 3}')
    text.AddText('#bf{Pb-Pb, #sqrt{s_{NN}} = 5.36 TeV}')
    text.AddText('#bf{#bar{p}-^{3}#bar{He}}')
    text.AddText(f'#bf{{#it{{N}}_{{raw}}(^{{4}}#bar{{Li}}) = {h_same_event_norm.Integral():.2f}}}')
    text.Draw()

    canvas.SaveAs('output/yield.pdf')

def compute_yield_error_me():

    h_same_event_norm = load_hist(
        HistLoadInfo('/home/galucia/antiLithium4/femto/output/AntiLithium4FitCFCATS.root',
        'dir/h_same_event_norm1.0')
    )

    outfile = TFile('output/yield_error.root', 'RECREATE')

    h_same_event_upper = h_same_event_norm.Clone('h_same_event_upper')
    h_same_event_lower = h_same_event_norm.Clone('h_same_event_lower')

    for i in range(1, h_same_event_norm.GetNbinsX() + 1):
        value = h_same_event_norm.GetBinContent(i)
        error = h_same_event_norm.GetBinError(i)

        h_same_event_upper.SetBinContent(i, value + error)
        h_same_event_lower.SetBinContent(i, value - error)

        
    integral_upper = h_same_event_upper.Integral()
    integral_lower = h_same_event_lower.Integral()
    yield_error = (integral_upper - integral_lower) / 2.0

    print(f'Yield error: {yield_error}')

    outfile.cd()
    h_same_event_norm.Write()
    h_same_event_upper.Write()
    h_same_event_lower.Write()

    yield_error_vector = TVectorD(1)
    yield_error_vector[0] = yield_error
    yield_error_vector.Write('yield_error')

    outfile.Close()



if __name__ == '__main__':
    draw_same_event_li()
    compute_yield_error_me()
    print('Done')