'''
    Classes for invariant mass studies
'''
import numpy as np

from ROOT import TH1F, TCanvas, TLine, TBox, TLegend, TF1
from ROOT import kGray, kOrange, kRed, gStyle

from .studies import StandaloneStudy

from torchic import HistLoadInfo
from torchic.core.histogram import load_hist
from torchic.utils.terminal_colors import TerminalColors as tc

import sys
sys.path.append('../..')
from framework.utils.root_setter import obj_setter

class InvariantMassStudy(StandaloneStudy):

    def __init__(self, config, outputFile, sameEvent=None, mixedEvent=None, **kwargs):
        '''
            Study to investigate the invariant mass distribution with different cuts.
        '''
        super().__init__(config, outputFile)
        self.opt = kwargs.get('opt', '')
        self.dir = self.outFile.mkdir('InvariantMass'+self.opt)

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

        self.hCorrection = None
        self.hSubtracted = None
        self.hPull = None
        self.hRatio = None

    def clone_same_event(self, sameEvent:TH1F) -> None:
        self.hSameEvent = sameEvent.Clone('hSame'+self.opt+'_invMass')

    def load_same_event(self, sameEventInfo:HistLoadInfo) -> None:
        self.hSameEvent = load_hist(sameEventInfo)
        self.hSameEvent.SetName('hSame'+self.opt+'_invMass')

    def clone_mixed_event(self, mixedEvent:TH1F) -> None:
        self.hMixedEvent = mixedEvent.Clone('hMixed'+self.opt+'_invMass')

    def load_mixed_event(self, mixedEventInfo:HistLoadInfo) -> None:
        self.hMixedEvent = load_hist(mixedEventInfo)
        self.hMixedEvent.SetName('hMixed'+self.opt+'_invMass')

    def set_same_event(self, sameEvent) -> None:
        if str(type(sameEvent)) == "<class 'ROOT.TH1F'>":                               self.clone_same_event(sameEvent)
        elif 'HistLoadInfo' in str(type(sameEvent)):                                    self.load_same_event(sameEvent)
        else:                                                                           raise ValueError('Type not supported')
    
    def set_mixed_event(self, mixedEvent) -> None:
        if str(type(mixedEvent)) == "<class 'ROOT.TH1F'>":                              self.clone_mixed_event(mixedEvent)
        elif 'HistLoadInfo' in str(type(mixedEvent)):                                   self.load_mixed_event(mixedEvent)
        else:                                                                           raise ValueError('Type not supported')

    def self_normalize(self) -> None:
        '''
            Normalize the sameEvent and the mixedEvent histograms.
        '''

        if self.hSameEvent:        
            low_edge = self.hSameEvent.GetBinLowEdge(1)
            high_edge = self.hSameEvent.GetBinLowEdge(self.hSameEvent.GetNbinsX()+1)
            sameEventIntegral = self.hSameEvent.Integral(1, self.hSameEvent.GetNbinsX(), 'width')
            self.hSameEvent.Scale((high_edge-low_edge)/sameEventIntegral, 'width')
        else:
            print(tc.GREEN+'[INFO]: '+tc.RESET+'No same event histogram provided')
        if self.hMixedEvent:
            low_edge = self.hMixedEvent.GetBinLowEdge(1)
            high_edge = self.hMixedEvent.GetBinLowEdge(self.hMixedEvent.GetNbinsX()+1)
            mixedEventIntegral = self.hMixedEvent.Integral(1, self.hMixedEvent.GetNbinsX(), 'width')
            self.hMixedEvent.Scale((high_edge-low_edge)/mixedEventIntegral, 'width')
        else:
            print(tc.GREEN+'[INFO]: '+tc.RESET+'No mixed event histogram provided')

    def normalize(self, low=3.78, high=3.89) -> None:
        '''
            Normalize the sameEvent and the mixedEvent histograms.
        '''

        if self.hMixedEvent and self.hSameEvent:
            low_bin = self.hMixedEvent.FindBin(low)
            low_edge = self.hMixedEvent.GetBinLowEdge(low_bin)
            high_bin = self.hMixedEvent.FindBin(high)
            high_edge = self.hMixedEvent.GetBinLowEdge(high_bin+1)
            sameEventIntegral = self.hSameEvent.Integral(low_bin, high_bin, 'width')
            mixedEventIntegral = self.hMixedEvent.Integral(low_bin, high_bin, 'width')
            self.hMixedEvent.Scale(sameEventIntegral/mixedEventIntegral)
        else:
            print(tc.GREEN+'[INFO]: '+tc.RESET+'No histogram provided')

    def correct_mixed_event(self, hCorrection: TH1F) -> None:
        '''
            Correct the mixed event histogram with a correction histogram.
            The correction histogram takes into account the effect of the interaction (obtained from the correlation function).
        '''
        if not self.hMixedEvent:
            print(tc.RED+'[ERROR]: '+tc.RESET+'No mixed event histogram provided')
            return
        
        self.hCorrection = hCorrection.Clone('hCorrection'+self.opt+'_invMass')

        for ibin in range(1, self.hMixedEvent.GetNbinsX()+1):
            mixedValue = self.hMixedEvent.GetBinContent(ibin)
            correctionValue = self.hCorrection.GetBinContent(ibin)
            self.hMixedEvent.SetBinContent(ibin, mixedValue*correctionValue)
            self.hMixedEvent.SetBinError(ibin, np.sqrt(mixedValue*correctionValue)) 

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
            tmp_hist = TH1F('hSame'+self.opt+'_invMass_tmp', 'hSame'+self.opt+'_invMass_tmp', len(bin_edges)-1, bin_edges)
            for ibin in range(1, self.hSameEvent.GetNbinsX()):
                tmp_hist.Fill(self.hSameEvent.GetBinCenter(ibin), self.hSameEvent.GetBinContent(ibin))
            for ibin in range(1, tmp_hist.GetNbinsX()+1):
                tmp_hist.SetBinError(ibin, np.sqrt(tmp_hist.GetBinContent(ibin)))
            self.hSameEvent = tmp_hist.Clone('hSame'+self.opt+'_invMass')
            del tmp_hist
            
        if self.hMixedEvent:
            tmp_hist = TH1F('hMixed'+self.opt+'_invMass_tmp', 'hMixed'+self.opt+'_invMass_tmp', len(bin_edges)-1, bin_edges)
            for ibin in range(1, self.hMixedEvent.GetNbinsX()):
                tmp_hist.Fill(self.hMixedEvent.GetBinCenter(ibin), self.hMixedEvent.GetBinContent(ibin))
            for ibin in range(1, tmp_hist.GetNbinsX()+1):
                tmp_hist.SetBinError(ibin, np.sqrt(tmp_hist.GetBinContent(ibin)))
            self.hMixedEvent = tmp_hist.Clone('hMixed'+self.opt+'_invMass')
            del tmp_hist

    # TODO: implement cut variation methods
        
    def bkg_subtraction(self) -> None:
        '''
            Subtract the background from the signal.
        '''
        if not self.hSameEvent or not self.hMixedEvent:
            print(tc.RED+'[ERROR]: '+tc.RESET+'No histograms provided')
            return

        self.hSubtracted = self.hSameEvent.Clone()
        self.hSubtracted.SetName('hSubtracted'+self.opt+'_invMass')
        self.hSubtracted.Reset()

        for bin in range(1, self.hSameEvent.GetNbinsX()+1):
            sameValue = self.hSameEvent.GetBinContent(bin)
            mixedValue = self.hMixedEvent.GetBinContent(bin)
            sameError = self.hSameEvent.GetBinError(bin)
            mixedError = self.hMixedEvent.GetBinError(bin)
            self.hSubtracted.SetBinContent(bin, sameValue-mixedValue)
            self.hSubtracted.SetBinError(bin, np.sqrt(sameError**2 + mixedError**2))

    def pull_distribution(self) -> None:
        '''
            Calculate the Pull distribution.
        '''
        if not self.hSubtracted:
            print(tc.RED+'[ERROR]: '+tc.RESET+'No subtracted histogram provided')
            return

        pol0 = TF1('pol0', 'pol0')
        pol0.SetParameter(0, 0)
        self.hSubtracted.Fit(pol0, 'MN+')

        param = pol0.GetParameter(0)
        error = pol0.GetParError(0)

        self.hPull = self.hSubtracted.Clone()
        self.hPull.SetName('hPull'+self.opt+'_invMass')
        self.hPull.Reset()

        for ibin in range(1, self.hSubtracted.GetNbinsX()+1):
            subtract_value = self.hSubtracted.GetBinContent(ibin)
            subtract_error = self.hSubtracted.GetBinError(ibin)
            pull = (subtract_value-param)/np.sqrt(subtract_error**2 + error**2)
            self.hPull.SetBinContent(ibin, pull)
            self.hPull.SetBinError(ibin, 0)

    def ratio_distribution(self) -> None:
        '''
            Calculate the ratio distribution.
        '''
        if not self.hSameEvent or not self.hMixedEvent:
            print(tc.RED+'[ERROR]: '+tc.RESET+'No histograms provided')
            return

        self.hRatio = self.hSameEvent.Clone()
        self.hRatio.SetName('hRatio'+self.opt+'_invMass')
        self.hRatio.Reset()

        for ibin in range(1, self.hSameEvent.GetNbinsX()+1):
            sameValue = self.hSameEvent.GetBinContent(ibin)
            mixedValue = self.hMixedEvent.GetBinContent(ibin)
            ratio = sameValue/mixedValue if mixedValue != 0 else 0
            ratioError = ratio*np.sqrt((self.hSameEvent.GetBinError(ibin)/sameValue)**2 + (self.hMixedEvent.GetBinError(ibin)/mixedValue)**2) if mixedValue != 0 and sameValue != 0 else 0
            self.hRatio.SetBinContent(ibin, ratio)
            self.hRatio.SetBinError(ibin, ratioError)

    def save(self, suffix='') -> None:
        '''
            Save the histograms in the output file.
        '''
        self.dir.cd()
        if self.hSameEvent:     self.hSameEvent.Write(self.hSameEvent.GetName()+suffix)
        if self.hMixedEvent:    self.hMixedEvent.Write(self.hMixedEvent.GetName()+suffix)
        if self.hCorrection:    self.hCorrection.Write(self.hCorrection.GetName()+suffix)
        if self.hSubtracted:    self.hSubtracted.Write(self.hSubtracted.GetName()+suffix)
        if self.hPull:          self.hPull.Write(self.hPull.GetName()+suffix)
        if self.hRatio:         self.hRatio.Write(self.hRatio.GetName()+suffix)
        canvas = TCanvas('canvas'+self.opt, 'canvas', 800, 600)
        canvas.cd()
        
        
        if self.hMixedEvent:         
            self.hMixedEvent.SetMarkerColor(kGray+2)
            self.hMixedEvent.Draw('same hist')
        if self.hSameEvent:     
            self.hSameEvent.SetMarkerColor(kOrange-3)
            self.hSameEvent.Draw('same')

        self.dir.cd()
        canvas.Write('canvas'+suffix)

    def produce_plot(self, outFilePath:str) -> None:
        '''
            Produce a plot with the invariant mass distribution.
        '''
        
        gStyle.SetOptStat(0)
        massLi = 3.74976
        widthLi = 0.006

        canvas = TCanvas('cInvMass'+self.opt, 'cInvMass', 800, 600)
        canvas.cd()
        obj_setter(self.hSubtracted, title='Invariant mass distribution;Invariant mass (GeV/#it{c}^{2});Counts', marker_style=20, marker_color=797)
        lineLi = TLine(massLi, -70, massLi, 70)
        obj_setter(lineLi, line_color=632, line_style=2)
        bandLi = TBox(massLi-widthLi/2, -70, massLi+widthLi/2, 70)
        obj_setter(bandLi, fill_color=632, fill_style=3004, fill_alpha=0.5)

        self.hSubtracted.Draw('e1 same')
        lineLi.Draw('same')
        bandLi.Draw('same')

        legend = TLegend(0.45, 0.12, 0.75, 0.3)
        legend.AddEntry(self.hSubtracted, 'data', 'p')
        legend.AddEntry(lineLi, '^{4}Li', 'l')
        legend.SetBorderSize(0)
        legend.SetTextSize(0.04)
        legend.Draw('same')

        self.dir.cd()
        canvas.Write()
        canvas.SaveAs(outFilePath)