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
        '_0.0',
        '_10.0',
        '_30.0',
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

def draw_invariant_mass_pad(h_same_event: TH1F, h_mixed_event: TH1F, pad_invmass: TPad, centrality: Centrality) -> None:

    pad_invmass.cd()
    pad_frame = pad_invmass.DrawFrame(
        h_same_event.GetXaxis().GetXmin(), h_same_event.GetMinimum() * 0.1,
        3.851, h_same_event.GetMaximum() * 1.2,
        '; m (p-^{3}He) (GeV/c^{2}); Counts (a.u.)'
    )
    
    pad_invmass.SetBottomMargin(0.005)
    pad_invmass.Draw()
    pad_invmass.SetLogy()
    pad_invmass.cd()
    
    h_same_event.SetTitle(';m (p-^{3}He) (GeV/c^{2}); Counts (a.u.)')
    h_same_event.SetLineColor(kOrange-3)
    h_same_event.SetLineWidth(2)
    h_same_event.SetMarkerStyle(22)
    h_same_event.SetMarkerColor(kOrange-3)
    h_same_event.SetMarkerSize(2)
    
    h_mixed_event.SetLineColor(kBlue+2)
    h_mixed_event.SetLineWidth(2)
    h_mixed_event.SetMarkerColor(kBlue+2)
    h_mixed_event.SetMarkerStyle(23)
    h_mixed_event.SetMarkerSize(2)

    li4_mass = 3.751 # GeV/c^2
    li4_width = 0.006 # GeV/c^2
    box = TBox(li4_mass - li4_width/2, 0, li4_mass + li4_width/2, h_same_event.GetMaximum() * 1.2)
    box.SetFillColorAlpha(kGreen+2, 0.2)
    
    h_same_event.Draw('e1 hist same')
    h_mixed_event.Draw('e1 hist same')
    box.Draw('same')

    legend  = TLegend(0.55, 0.05, 0.89, 0.3)
    legend.SetBorderSize(0)
    legend.SetFillColor(0)
    legend.AddEntry(h_same_event, 'Same event', 'lep')
    legend.AddEntry(h_mixed_event, 'Mixed event', 'lep')
    legend.AddEntry(box, '^{4}Li 1#sigma band', 'f')
    legend.Draw('same')

    watermark = get_alice_watermark(0.55, 0.35, 0.89, 0.6)
    watermark.Draw('same')
    centrality_text = TPaveText(0.55, 0.3, 0.89, 0.35, 'NDC')
    centrality_text.SetBorderSize(0)
    centrality_text.SetFillColor(0)
    centrality_text.SetTextFont(42)
    text = 'Centrality 0 - 50 %' if centrality == Centrality.INTEGRATED else \
        'Centrality: ' + SUFFIX_DICT['undot'][centrality.value].removeprefix('_cent').replace('_', ' - ') + '%'
    centrality_text.AddText(text)
    centrality_text.Draw('same')

    return box, watermark, legend, centrality_text


def draw_pull_pad(pull: TH1F, pad_pull: TPad) -> None:
    '''
        Draw the pull pad with the zero line.
    '''
    pad_pull.SetTopMargin(0.0)
    pad_pull.SetBottomMargin(0.25)
    pad_pull.cd()
    pull_frame = pad_pull.DrawFrame(
        pull.GetXaxis().GetXmin(), -5, 
        3.851, 5,
        ';m (p-^{3}He) (GeV/c^{2});n#sigma'
    )


    pull_frame.GetXaxis().SetLabelSize(.08)
    pull_frame.GetYaxis().SetLabelSize(.08)
    pull_frame.GetXaxis().SetTitleSize(.08)
    pull_frame.GetYaxis().SetTitleSize(.1)
    pull_frame.GetYaxis().SetTitleOffset(0.5)
    pull.SetMarkerStyle(20)
    pull.SetMarkerSize(2)
    pull.SetMarkerColor(kGreen+2)
    
    zero_line = TF1('zero_line', '0', pull.GetXaxis().GetXmin(), pull.GetXaxis().GetXmax())
    zero_line.SetLineColor(kGray+1)
    zero_line.SetLineStyle(2)
    zero_line.SetLineWidth(2)

    li4_mass = 3.751 # GeV/c^2
    li4_width = 0.006 # GeV/c^2
    box = TBox(li4_mass - li4_width/2, -5, li4_mass + li4_width/2, 5)
    pull.GetYaxis().Print()
    print(f'Pull Y axis range: {pull.GetYaxis().GetXmin()} - {pull.GetYaxis().GetXmax()}')
    box.SetFillColorAlpha(kGreen+2, 0.2)
    
    pad_pull.cd()
    pull.Draw('E1 same')
    zero_line.Draw('same')
    box.Draw('same')

    return zero_line, box

def draw_canvas(h_same_event, h_mixed_event, pull, pdf_outpath, centrality: Centrality) -> None:

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
    
    _box, _watermark, _legend, _centrality_text = \
        draw_invariant_mass_pad(h_same_event, h_mixed_event, pad_invmass, centrality)
    _zero_line, _pull_box = draw_pull_pad(pull, pad_pull)
    
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

def load_hists(input_file, args, centrality: Centrality):
    
    file = TFile.Open(input_file)
    sign = args.sign
    suffix = SUFFIX_DICT['dot'][centrality.value]

    h_same_event = file.Get(f'InvariantMass{sign}/hSame_invMass{suffix}')
    h_same_event.SetDirectory(0)  # Detach from file to avoid issues with deletion
    h_mixed_event = file.Get(f'InvariantMass{sign}/hMixed_invMass{suffix}')
    h_mixed_event.SetDirectory(0)  # Detach from file to avoid issues with deletion

    # renormalise
    integral_same_event = h_same_event.Integral()
    integral_mixed_event = h_mixed_event.Integral()
    h_mixed_event.Scale(integral_same_event/integral_mixed_event)

    if centrality != Centrality.INTEGRATED:
        h_same_event.Rebin()
        h_mixed_event.Rebin()

    h_same_event.GetXaxis().SetRangeUser(h_same_event.GetXaxis().GetXmin(), 3.851)
    h_mixed_event.GetXaxis().SetRangeUser(h_mixed_event.GetXaxis().GetXmin(), 3.851)

    pull = perform_gaussian_pull(h_same_event, h_mixed_event)
    pull.SetDirectory(0)  # Detach from file to avoid issues with deletion

    file.Close()
    return h_same_event, h_mixed_event, pull

def main(args: argparse.Namespace):

    #watermark = get_alice_watermark(0.6, 0.7, 0.95, 0.9)
    
    pdf_outpath = 'output/invmass_comparison.pdf'
    blank_canvas = TCanvas('blank_canvas', 'blank_canvas', 1000, 1000)
    blank_canvas.Print(pdf_outpath + '(')  # Create a blank PDF file
    
    for centrality in Centrality:
        print(f'Processing centrality: {centrality.name}')
        h_same_event, h_mixed_event, pull = load_hists('/home/galucia/antiLithium4/analysis/output/PbPb/studies_us.root', args, centrality)
        draw_canvas(h_same_event, h_mixed_event, pull, pdf_outpath, centrality)
    
    blank_canvas.Print(pdf_outpath + ')')  # Close the PDF file

    
    # Save the watermark
    #watermark.Draw()
    
    print(f'Figures saved to {pdf_outpath}')


if __name__ == '__main__':
    args = parse_args()
    main(args)
