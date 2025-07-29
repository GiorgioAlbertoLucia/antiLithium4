
from ROOT import TFile, TCanvas, TIter, kRed, kGreen, TLegend, TPaveText

infile = TFile.Open('OutputTreeHe3NEW_noLB25MC.root')

canvas = infile.Get('DCAsamePt_MC_Data/DCA_SecvsPrim_2.35')
primitives = canvas.GetListOfPrimitives()

h_prim, h_sec = None, None

next_item = TIter(primitives)
while (item := next_item()):
    if item.GetName() == 'hist_histDCAxyPrim_2.35':
        h_prim = item
        print('Found histogram for primaries')
    elif item.GetName() == 'hist_histDCAxySec_2.35':
        h_sec = item
        print('Found histogram for secondaries')

canvas = TCanvas('canvas', 'canvas', 800, 600)
h_prim.SetLineColor(kRed)
h_prim.SetLineWidth(2)

h_sec.SetLineColor(kGreen+2)
h_sec.SetLineWidth(2)


h_prim.Draw('hist')
h_sec.Draw('hist same')

legend = TLegend(0.6, 0.65, 0.88, 0.85)
legend.SetBorderSize(0)
legend.SetFillColor(0)
legend.AddEntry(h_prim, 'Primaries', 'l')
legend.AddEntry(h_sec, 'Secondaries', 'l')
legend.Draw()

text = TPaveText(0.15, 0.55, 0.4, 0.85, 'NDC')
text.SetFillColor(0)
text.SetBorderSize(0)
text.AddText('WIP in the ALICE')
text.AddText('group in Turin')
text.AddText('#bf{ALICE Run 3}')
text.AddText('#bf{pp #sqrt{s} = 13.6 TeV}')
text.AddText('#bf{2.1 < #it{p}_{T} < 2.6 GeV/#it{c}}')
canvas.cd()
text.Draw()

canvas.SaveAs('DCAxyPrimSec_2.35.pdf')

