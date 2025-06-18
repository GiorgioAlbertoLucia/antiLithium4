'''
    Simple script to prepare invariant mass comparison for same event and mixed event. The pull is also produced.
'''

import argparse
from enum import Enum
import numpy as np

from ROOT import TFile, TH1F, TCanvas, gStyle, TPaveText, TPad, TF1, TBox, TLegend
from ROOT import kOrange, kBlue, kGreen, kGray

class Centrality(Enum):
    CENT_0_10 = 0
    CENT_10_30 = 1
    CENT_30_50 = 2
    INTEGRATED = 3

SUFFIX_DICT = {
    'dot': 
    [
        '_cent0.0_10.0',
        '_cent10.0_30.0',
        '_cent30.0_50.0',
        '',
    ],
    'undot':
    [
        '_cent0_10',
        '_cent10_30',
        '_cent30_50',
        '',
    ]
}

def get_alice_watermark(xmin, ymin, xmax, ymax):

    '''
        Get the Alice watermark
    '''
    watermark = TPaveText(xmin, ymin, xmax, ymax, 'NDC')
    watermark.SetBorderSize(0)
    watermark.SetFillColor(0)
    watermark.AddText('This work')
    watermark.AddText('#bf{ALICE, Run 3}')
    watermark.AddText('#bf{Pb#topbarPb, #sqrt{s_{NN}} = 5.36 TeV}')
    return watermark

def parse_args():
    parser = argparse.ArgumentParser(description='Prepare invariant mass comparison for same event and mixed event.')
    parser.add_argument('--sign', type=str, default='Anti', choices=['Matter', 'Anti'], help='Sign of the mother particle (Matter or Anti)')
    return parser.parse_args()

def draw_correlation_function_pad(h_correlation_experimental: TH1F, h_correlation_theory: TH1F, pad_invmass: TPad, centrality: Centrality) -> None:

    pad_invmass.cd()
    pad_frame = pad_invmass.DrawFrame(
        h_correlation_experimental.GetXaxis().GetXmin(), 0.4,
        0.8, h_correlation_experimental.GetMaximum() * 1.2,
        ';#it{k*} (GeV/#it{c});C(#it{k*})'
    )
    
    pad_invmass.SetBottomMargin(0.005)
    pad_invmass.Draw()
    pad_invmass.cd()
    
    h_correlation_experimental.SetTitle(';m (p-^{3}He) (GeV/c^{2}); Counts (a.u.)')
    h_correlation_experimental.SetLineColor(kOrange-3)
    h_correlation_experimental.SetLineWidth(2)
    h_correlation_experimental.SetMarkerStyle(22)
    h_correlation_experimental.SetMarkerColor(kOrange-3)
    h_correlation_experimental.SetMarkerSize(2)
    
    h_correlation_theory.SetLineColor(kBlue+2)
    h_correlation_theory.SetLineWidth(2)
    h_correlation_theory.SetMarkerColor(kBlue+2)
    h_correlation_theory.SetMarkerStyle(23)
    h_correlation_theory.SetMarkerSize(1)

    li4_peak = 0.081 # GeV/c
    li4_width = 0.014 * 3 # GeV/c
    #box = TBox(li4_peak - li4_width/2, 0, li4_peak + li4_width/2, h_correlation_experimental.GetMaximum() * 1.2)
    #box.SetFillColorAlpha(kGreen+2, 0.2)
    
    h_correlation_experimental.Draw('e1 same')
    h_correlation_theory.Draw('c same')
    #box.Draw('same')

    legend  = TLegend(0.55, 0.45, 0.89, 0.6)
    legend.SetBorderSize(0)
    legend.SetFillColor(0)
    legend.AddEntry(h_correlation_experimental, 'Data', 'lep')
    legend.AddEntry(h_correlation_theory, 'CATS', 'lep')
    #legend.AddEntry(box, '^{4}Li 3#sigma band', 'f')
    legend.Draw('same')

    watermark = get_alice_watermark(0.55, 0.65, 0.89, 0.89)
    watermark.Draw('same')
    centrality_text = TPaveText(0.55, 0.6, 0.89, 0.65, 'NDC')
    centrality_text.SetBorderSize(0)
    centrality_text.SetFillColor(0)
    centrality_text.SetTextFont(42)
    text = 'Centrality 0 - 50 %' if centrality == Centrality.INTEGRATED else \
        'Centrality: ' + SUFFIX_DICT['undot'][centrality.value].removeprefix('_cent').replace('_', ' - ') + '%'
    centrality_text.AddText(text)
    centrality_text.Draw('same')

    return watermark, legend, centrality_text

def draw_pull_pad(pull: TH1F, pad_pull: TPad) -> None:
    '''
        Draw the pull pad with the zero line.
    '''
    pad_pull.SetTopMargin(0.0)
    pad_pull.SetBottomMargin(0.25)
    pad_pull.cd()
    pull_frame = pad_pull.DrawFrame(
        pull.GetXaxis().GetXmin(), -5, 
        0.8, 5,
        ';#it{k*} (GeV/#it{c});n#sigma'
    )

    pull_frame.GetXaxis().SetLabelSize(.08)
    pull_frame.GetYaxis().SetLabelSize(.08)
    pull_frame.GetXaxis().SetTitleSize(.08)
    pull_frame.GetYaxis().SetTitleSize(.1)
    pull_frame.GetYaxis().SetTitleOffset(0.5)
    pull.SetMarkerStyle(20)
    pull.SetMarkerSize(2)
    pull.SetMarkerColor(kGreen+2)

    for ibin in range(1, pull.GetNbinsX()+1):
        print(f'Bin {ibin}: {pull.GetBinContent(ibin)} ± {pull.GetBinError(ibin)}')
    
    zero_line = TF1('zero_line', '0', pull.GetXaxis().GetXmin(), pull.GetXaxis().GetXmax())
    zero_line.SetLineColor(kGray+1)
    zero_line.SetLineStyle(2)
    zero_line.SetLineWidth(2)

    li4_peak = 0.081 # GeV/c
    li4_width = 0.014 * 3 # GeV/c
    #box = TBox(li4_peak - li4_width/2, -5, li4_peak + li4_width/2, 5)
    #box.SetFillColorAlpha(kGreen+2, 0.2)
    
    pad_pull.cd()
    pull.Draw('E1 same')
    zero_line.Draw('same')
    #box.Draw('same')

    return zero_line#, box

def draw_canvas(h_correlation_experimental, h_correlation_theory, pull, pdf_outpath, centrality: Centrality) -> None:

    '''
        Draw the canvas with invariant mass and pull.
    '''

    gStyle.SetOptStat(0)
    
    canvas = TCanvas('cf_canvas', 'cf_canvas', 1000, 1000)
    pad_invmass = TPad('pad_invmass', 'pad_invmass', 0, 0.3, 1, 1)
    pad_pull = TPad('pad_pull', 'pad_pull', 0, 0, 1, 0.3)
    
    canvas.cd()
    pad_invmass.Draw()
    pad_pull.Draw()
    
    _watermark, _legend, _centrality_text = \
        draw_correlation_function_pad(h_correlation_experimental, h_correlation_theory, pad_invmass, centrality)
    _zero_line = draw_pull_pad(pull, pad_pull)
    
    canvas.SetTopMargin(0.05)
    canvas.SetBottomMargin(0.15)
    canvas.Print(pdf_outpath)

def perform_gaussian_pull(h1, h2) -> TH1F:

    pull = h1.Clone('pull')
    pull.Reset()

    for ibin in range(1, h1.GetNbinsX() + 1):
        
        error = np.sqrt(h2.GetBinError(ibin)**2 + h1.GetBinError(ibin)**2)
        error = error if error > 0 else 1e-6
        pull_value = (h1.GetBinContent(ibin) - h2.GetBinContent(ibin)) / error
        pull.SetBinContent(ibin, pull_value)
        pull.SetBinError(ibin, 1.)
    
    pull.SetTitle(';m (p-^{3}He) (GeV/c^{2});n#sigma')
    return pull

def load_hists(input_file_cf, input_file_theory, args, centrality: Centrality):
    
    sign = args.sign
    suffix = SUFFIX_DICT['dot'][centrality.value]

    file = TFile.Open(input_file_cf)
    h_correlation_experimental = file.Get(f'Correlation{sign}/hCorrelation_kstar{suffix}')
    h_correlation_experimental.SetDirectory(0)  # Detach from file to avoid issues with deletion
    file.Close()
    
    file = TFile.Open(input_file_theory)
    h_correlation_theory = file.Get(f'hHe3_p_Coul_CF_US')
    h_correlation_theory.SetDirectory(0)  # Detach from file to avoid issues with deletion
    h_correlation_theory.GetXaxis().SetRangeUser(h_correlation_experimental.GetXaxis().GetXmin(), h_correlation_experimental.GetXaxis().GetXmax())
    file.Close()

    h_correlation_experimental.Rebin()
    h_correlation_experimental.Scale(1/2 )

    h_correlation_theory_for_pull = h_correlation_theory.Clone('h_correlation_theory_for_pull')
    h_correlation_theory_for_pull.Rebin(20)
    h_correlation_theory_for_pull.Scale(1/20)

    pull = perform_gaussian_pull(h_correlation_experimental, h_correlation_theory_for_pull)
    pull.SetDirectory(0)  # Detach from file to avoid issues with deletion

    return h_correlation_experimental, h_correlation_theory, pull



def main(args: argparse.Namespace):

    #watermark = get_alice_watermark(0.6, 0.7, 0.95, 0.9)
    
    pdf_outpath = 'output/correlation_us.pdf'
    blank_canvas = TCanvas('blank_canvas', 'blank_canvas', 1000, 1000)
    blank_canvas.Print(pdf_outpath + '(')  # Create a blank PDF file

    infile_experimental = '/home/galucia/antiLithium4/analysis/output/PbPb/studies_us.root'
    
    for centrality in Centrality:
        print(f'Processing centrality: {centrality.name}')
        suffix = SUFFIX_DICT['undot'][centrality.value]
        infile_theory = f'/home/galucia/antiLithium4/analysis/output/CATS/CATS{suffix}_new.root'
        h_correlation_experimental, h_correlation_theory, pull = load_hists(
            infile_experimental, infile_theory, args, centrality)
        draw_canvas(h_correlation_experimental, h_correlation_theory, pull, pdf_outpath, centrality)
    
    blank_canvas.Print(pdf_outpath + ')')  # Close the PDF file

    
    # Save the watermark
    #watermark.Draw()
    
    print(f'Figures saved to {pdf_outpath}')


if __name__ == '__main__':
    args = parse_args()
    main(args)
