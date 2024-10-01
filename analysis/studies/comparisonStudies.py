'''
    Classes for correlation function studies
'''
import numpy as np
from ROOT import TH1F, TCanvas, TLine, TLegend
from ROOT import kOrange, kGray

from .studies import StandaloneStudy

import sys

sys.path.append('../..')
from framework.src.hist_info import HistLoadInfo
from framework.src.hist_handler import HistHandler
from framework.utils.terminal_colors import TerminalColors as tc
from framework.utils.root_setter import obj_setter

class ComparisonStudy(StandaloneStudy):

    def __init__(self, config, sameEvent=None, mixedEvent=None):
        '''
            Study to investigate the invariant mass distribution with different cuts.
        '''
        super().__init__(config)
        self.dir = ComparisonStudy.outFile_shared.mkdir('Comparison')

        if sameEvent:
            self.set_same_event(sameEvent)
        else:
            print(tc.MAGENTA+'[WARNING]: '+tc.RESET+'No same event histogram provided')
            self.hSameEvent = None
        if mixedEvent:
            self.set_mixed_event(mixedEvent)
        else:
            print(tc.MAGENTA+'[WARNING]: '+tc.RESET+'No mixedEvent histogram provided')
            self.hMixedEvent = None

        self.hCorrelation = None

    def clone_same_event(self, sameEvent:TH1F) -> None:
        self.hSameEvent = sameEvent.Clone('hSame_'+sameEvent.GetName())

    def load_same_event(self, sameEventInfo:HistLoadInfo) -> None:
        self.hSameEvent = HistHandler.loadHist(sameEventInfo)
        self.hSameEvent.SetName('hSame_'+self.hSameEvent.GetName())

    def clone_mixed_event(self, mixedEvent:TH1F) -> None:
        self.hMixedEvent = mixedEvent.Clone('hMixed_'+mixedEvent.GetName())

    def load_mixed_event(self, mixedEventInfo:HistLoadInfo) -> None:
        self.hMixedEvent = HistHandler.loadHist(mixedEventInfo)
        self.hMixedEvent.SetName('hMixed_'+self.hMixedEvent.GetName())

    def set_same_event(self, sameEvent) -> None:
        if str(type(sameEvent)) == "<class 'ROOT.TH1F'>":                               self.clone_same_event(sameEvent)
        elif str(type(sameEvent)) == "<class 'framework.src.hist_info.HistLoadInfo'>":  self.load_same_event(sameEvent)
        else:                                                                           raise ValueError('Type not supported')
    
    def set_mixed_event(self, mixedEvent) -> None:
        if str(type(mixedEvent)) == "<class 'ROOT.TH1F'>":                              self.clone_mixed_event(mixedEvent)
        elif str(type(mixedEvent)) == "<class 'framework.src.hist_info.HistLoadInfo'>": self.load_mixed_event(mixedEvent)
        else:                                                                           raise ValueError('Type not supported')

    def self_normalize(self) -> None:
        '''
            Normalize the sameEvent and the mixedEvent histograms.
        '''

        if self.hSameEvent:        
            low_edge = self.hSameEvent.GetBinLowEdge(1)
            high_edge = self.hSameEvent.GetBinLowEdge(self.hSameEvent.GetNbinsX()+1)
            sameEventIntegral = self.hSameEvent.Integral(1, self.hSameEvent.GetNbinsX(), 'width')
            self.hSameEvent.Scale((high_edge-low_edge)/sameEventIntegral)
            for ibin in range(1, self.hSameEvent.GetNbinsX()):
                ivalue = self.hSameEvent.GetBinContent(ibin)
                ierror = np.sqrt((high_edge-low_edge)*ivalue/sameEventIntegral + ivalue*ivalue/sameEventIntegral)
                self.hSameEvent.SetBinError(ibin, ierror)
        else:
            print(tc.GREEN+'[INFO]: '+tc.RESET+'No same event histogram provided')
        if self.hMixedEvent:
            low_edge = self.hMixedEvent.GetBinLowEdge(1)
            high_edge = self.hMixedEvent.GetBinLowEdge(self.hMixedEvent.GetNbinsX()+1)
            mixedEventIntegral = self.hMixedEvent.Integral(1, self.hMixedEvent.GetNbinsX(), 'width')
            self.hMixedEvent.Scale((high_edge-low_edge)/mixedEventIntegral)
            for ibin in range(1, self.hMixedEvent.GetNbinsX()):
                ivalue = self.hMixedEvent.GetBinContent(ibin)
                ierror = np.sqrt((high_edge-low_edge)*ivalue/mixedEventIntegral + ivalue*ivalue/mixedEventIntegral)
                self.hMixedEvent.SetBinError(ibin, ierror)
        else:
            print(tc.GREEN+'[INFO]: '+tc.RESET+'No mixed event histogram provided')


    def scale_mixing(self, scale_factor=1.0) -> None:
        '''
            Normalize the sameEvent and the mixedEvent histograms.
        '''

        if self.hMixedEvent:
            self.hMixedEvent.Scale(scale_factor)
        else:
            print(tc.GREEN+'[INFO]: '+tc.RESET+'No histogram provided')

    def rebin(self, rebin_factor:int=2) -> None:
        '''
            Rebin the sameEvent and the mixedEvent histograms.
        '''
        if self.hSameEvent:     self.hSameEvent.Rebin(rebin_factor)
        if self.hMixedEvent:    self.hMixedEvent.Rebin(rebin_factor)

    def custom_binning(self, bin_edges:np.ndarray) -> None:
        '''
            Define a custom binning for the histograms
        '''
        if self.hSameEvent:     
            tmp_hist = TH1F('hSame_tmp', 'hSame_tmp', len(bin_edges)-1, bin_edges)
            for ibin in range(1, self.hSameEvent.GetNbinsX()+1):
                tmp_hist.Fill(self.hSameEvent.GetBinCenter(ibin), self.hSameEvent.GetBinContent(ibin))
            for ibin in range(1, tmp_hist.GetNbinsX()):
                tmp_hist.SetBinError(ibin, np.sqrt(tmp_hist.GetBinContent(ibin)))
            self.hSameEvent = tmp_hist.Clone('hSame_'+self.hSameEvent.GetName())
            del tmp_hist
            
        if self.hMixedEvent:
            tmp_hist = TH1F('hMixed_tmp', 'hMixed_tmp', len(bin_edges)-1, bin_edges)
            for ibin in range(1, self.hMixedEvent.GetNbinsX()+1):
                tmp_hist.Fill(self.hMixedEvent.GetBinCenter(ibin), self.hMixedEvent.GetBinContent(ibin))
            for ibin in range(1, tmp_hist.GetNbinsX()):
                tmp_hist.SetBinError(ibin, np.sqrt(tmp_hist.GetBinContent(ibin)))
            self.hMixedEvent = tmp_hist.Clone('hMixed_'+self.hMixedEvent.GetName())
            del tmp_hist

    def save(self, suffix='') -> None:

        self.dir.cd()
        if self.hSameEvent:     self.hSameEvent.Write(self.hSameEvent.GetName()+suffix)
        if self.hMixedEvent:    self.hMixedEvent.Write(self.hMixedEvent.GetName()+suffix)
        canvas = TCanvas('canvas_comparison_'+self.hSameEvent.GetName()+'_'+suffix, 'canvas', 800, 600)
        canvas.cd()
        
        if self.hSameEvent:     
            obj_setter(self.hSameEvent, marker_style=20, marker_color=797)
            self.hSameEvent.Draw('e1 same')
        if self.hMixedEvent:         
            obj_setter(self.hMixedEvent, marker_style=20, marker_color=857)
            self.hMixedEvent.Draw('e1 same')

        legend = TLegend(0.15, 0.6, 0.49, 0.8)
        if self.hSameEvent:     legend.AddEntry(self.hSameEvent, 'Same event', 'p')
        if self.hMixedEvent:    legend.AddEntry(self.hMixedEvent, 'Mixed event', 'p')
        legend.SetBorderSize(0)
        legend.SetTextSize(0.04)
        legend.Draw('same')
        
        self.dir.cd()
        canvas.Write()

    def produce_plot(self, output_pdf:str) -> None:
        '''
            Produce the plots
        '''

        canvas = TCanvas('canvas', 'canvas', 800, 600)
        canvas.cd()
        #hframe = canvas.DrawFrame(0.0, 0., 1.0, 2., ';k* (GeV/#it{c}); C(k*)')
        obj_setter(self.hSameEvent, marker_style=20, marker_color=797)
        self.hSameEvent.Draw('e1 same')
        obj_setter(self.hMixedEvent, marker_style=20, marker_color=857)
        self.hMixedEvent.Draw('e1 same')
        
        legend = TLegend(0.15, 0.6, 0.49, 0.8)
        legend.AddEntry(self.hSameEvent, 'Same event', 'p')
        legend.AddEntry(self.hMixedEvent, 'Mixed event', 'p')
        legend.SetBorderSize(0)
        legend.SetTextSize(0.04)
        legend.Draw('same')
        
        canvas.SaveAs(output_pdf)
        print(tc.GREEN+'[INFO]: '+tc.RESET+'Saved plot to '+output_pdf)