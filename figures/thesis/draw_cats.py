
from ROOT import TGraph, TFile, kOrange, TCanvas, gStyle, kBlue, kGreen, TLegend, TPaveText

graphs = []
for file, name, color in zip(['/home/galucia/antiLithium4/analysis/output/CATS/CATS_cent0_10_converted.root',
                              '/home/galucia/antiLithium4/analysis/output/CATS/CATS_cent10_30_converted.root',
                              '/home/galucia/antiLithium4/analysis/output/CATS/CATS_cent30_50_converted.root',], 
                              ['0-10%', '10-30%', '30-50%'], [kOrange-3, kBlue-3, kGreen-3]):

    infile = TFile.Open(file)
    hist = infile.Get("hHe3_p_Coul_CF")

    graph = TGraph(hist.GetNbinsX())
    for ibin in range(1, hist.GetNbinsX() + 1):    
        graph.SetPoint(ibin - 1, hist.GetBinCenter(ibin), hist.GetBinContent(ibin))

    graph.SetName(name)
    graph.SetLineColor(color)
    graph.SetLineWidth(2)
    print(f"Graph {name} has {graph.GetN()} points")

    graphs.append(graph)

gStyle.SetOptStat(0)
canvas = TCanvas("canvas", "canvas", 800, 600)
canvas.SetLeftMargin(0.15)
frame = canvas.DrawFrame(0, 0.2, 0.4, 1.2, "; #it{k}* (GeV/#it{c}); C(#it{k}*)") 

text = TPaveText(0.4, 0.6, 0.5, 0.64, "NDC")
text.SetFillColor(0)
text.SetBorderSize(0)
text.SetTextSize(0.04)
text.AddText("This work")

legend = TLegend(0.4, 0.4, 0.8, 0.6, '#bf{CATS}')
legend.SetFillColor(0)
legend.SetBorderSize(0)
legend.SetTextSize(0.04)

canvas.cd()
for igraph, graph in enumerate(graphs):
    graph.Draw("L SAME")
    legend.AddEntry(graph, graph.GetName(), "l")

legend.Draw()
text.Draw()

canvas.SaveAs("CATS.pdf")