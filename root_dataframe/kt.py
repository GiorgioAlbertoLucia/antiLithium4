
import ROOT

ROOT.gStyle.SetOptStat(0)
canvas = ROOT.TCanvas("canvas", "canvas", 800, 600)
canvas.SetLeftMargin(0.15)
frame = canvas.DrawFrame(0.8, 0, 3, 27100, "; #it{m}_{T} (GeV/#it{c}^{2}); Counts")

infile = ROOT.TFile.Open("output/same_event.root")
kt010 = infile.Get("kstarAntimatter/hMt010Antimatter")
kt1030 = infile.Get("kstarAntimatter/hMt1030Antimatter")
kt3050 = infile.Get("kstarAntimatter/hMt3050Antimatter")

kt010.SetLineWidth(2)
kt1030.SetLineWidth(2)
kt3050.SetLineWidth(2)
kt010.SetLineColor(ROOT.kRed-2)
kt1030.SetLineColor(ROOT.kBlue-2)
kt3050.SetLineColor(ROOT.kGreen-2)

means = [kt010.GetMean(), kt1030.GetMean(), kt3050.GetMean()]
cents = ["0-10%", "10-30%", "30-50%"]

text = ROOT.TPaveText(0.6, 0.56, 0.8, 0.88, 'ndc')
text.SetFillColor(0)
text.SetBorderSize(0)
text.SetTextSize(0.03)
text.AddText('This work')
text.AddText("#bf{ALICE Run 3}")
text.AddText("#bf{Pb-Pb #sqrt{s_{NN}} = 5.36 TeV}")
for cent, mean in zip(cents, means):
    text.AddText(f"#bf{{#LT #it{{m}}_{{T}} #GT ({cent}) = {mean:.2f} GeV/#it{{c}}^{{2}}}}")

kt010.SetTitle("; #it{m}_{T} (GeV/#it{c}); Counts")
kt010.Draw("HIST same")
kt1030.Draw("HIST SAME")
kt3050.Draw("HIST SAME")

legend = ROOT.TLegend(0.64, 0.4, 0.8, 0.56)
legend.SetFillColor(0)
legend.SetBorderSize(0)
legend.AddEntry(kt010, "0-10%", "l")
legend.AddEntry(kt1030, "10-30%", "l")
legend.AddEntry(kt3050, "30-50%", "l")

legend.Draw()
text.Draw()
canvas.SaveAs("output/mt_histograms.pdf")

mean_050 = 0
total = 0
for hist in [kt010, kt1030, kt3050]:
    mean_050 += hist.GetMean() * hist.GetEntries()
    total += hist.GetEntries()

mean_050 /= total if total > 0 else 1
print(f"Mean kT for 0-50% centrality: {mean_050:.2f} GeV/c")