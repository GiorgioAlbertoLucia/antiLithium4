# code to convert a TGraph to a TH1F

from ROOT import TGraph, TH1F, TFile

def graph_to_hist(graph):
    hist = TH1F(f'h_{graph.GetName()}', graph.GetTitle(), graph.GetN(), graph.GetX()[0], graph.GetX()[graph.GetN()-1])
    for i in range(graph.GetN()):
        hist.SetBinContent(i+1, graph.GetY()[i])
    return hist

def scale_x_axis(hist, scale_factor):
    hist.SetTitle(hist.GetTitle() + f" ({hist.GetXaxis().GetTitle()} in GeV)")
    hist.GetXaxis().SetTitle(hist.GetXaxis().GetTitle() + " (GeV)")
    hist.GetXaxis().SetLimits(hist.GetXaxis().GetXmin() / scale_factor, hist.GetXaxis().GetXmax() / scale_factor)
    return hist

def main():

    infile_path = "/home/galucia/antiLithium4/analysis/output/fOutputSc.root"
    ingraph_name = "grHe3_p_CF_CoulAttr"
    outfile_path = "/home/galucia/antiLithium4/analysis/output/hist_fOutputSc.root"

    infile = TFile(infile_path, "READ")
    graph = infile.Get(ingraph_name)
    hist = graph_to_hist(graph)
    scale_x_axis(hist, scale_factor=1000)

    outfile = TFile(outfile_path, "RECREATE")
    graph.Write()
    hist.Write()
    outfile.Close()

def main2(infile_path: str, outfile: TFile): 

    infile = TFile(infile_path, "READ")
    hist = infile.Get("hHe3_p_Coul_CF")
    scale_x_axis(hist, scale_factor=1000)
    if 'LS' in infile_path:
        hist.SetName("hHe3_p_Coul_CF_LS")
    else:
        hist.SetName("hHe3_p_Coul_CF_US")
    outfile.cd()
    hist.Write()

def convert_and_rebin(infile_path:str, outfile: TFile):

    infile = TFile(infile_path, "READ")
    hist = infile.Get("hHe3_p_Coul_InvMass")
    scale_x_axis(hist, scale_factor=1000)
    hist.SetName("hHe3_p_Coul_InvMass")
    hist.Rebin(4)
    hist.Scale(1./4)

    mmin = 3.747
    mmax = 3.851    
    bin_width =  0.002
    nbins = int((mmax - mmin) / bin_width)
    tmp_hist = TH1F("tmp_hist", "tmp_hist", nbins, mmin, mmax)
    for i in range(1, nbins+1):
        tmp_hist.SetBinContent(i, hist.GetBinContent(i))
    outfile.cd()
    tmp_hist.Write("hHe3_p_Coul_InvMass")




if __name__ == "__main__":
    
    #main()

    outfile = TFile("/home/galucia/antiLithium4/analysis/output/CATS/CATS_scaled.root", "RECREATE")
    main2("/home/galucia/antiLithium4/analysis/output/CATS/CATS_CF_LS.root", outfile)
    main2("/home/galucia/antiLithium4/analysis/output/CATS/CATS_CF_US.root", outfile)
    convert_and_rebin("/home/galucia/antiLithium4/analysis/output/CATS/CATS_CF_LS.root", outfile)
    outfile.Close()