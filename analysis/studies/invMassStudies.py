'''
    Classes for invariant mass studies
'''
import numpy as np

from ROOT import TH1F

from .studies import StandaloneStudy

import sys
sys.path.append('../..')
from framework.src.hist_info import HistLoadInfo
from framework.src.hist_handler import HistHandler
from framework.utils.terminal_colors import TerminalColors as tc

class InvariantMassStudy(StandaloneStudy):

    def __init__(self, config, sameEvent=None, mixedEvent=None):
        '''
            Study to investigate the invariant mass distribution with different cuts.
        '''
        super().__init__(config)
        self.dir = InvariantMassStudy.outFile_shared.mkdir('InvariantMass')

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

    def clone_same_event(self, sameEvent:TH1F) -> None:
        self.hSameEvent = sameEvent.Clone('hSame_invMass')

    def load_same_event(self, sameEventInfo:HistLoadInfo) -> None:
        self.hSameEvent = HistHandler.loadHist(sameEventInfo)
        self.hSameEvent.SetName('hSame_invMass')

    def clone_mixed_event(self, mixedEvent:TH1F) -> None:
        self.hMixedEvent = mixedEvent.Clone('hMixed_invMass')

    def load_mixed_event(self, mixedEventInfo:HistLoadInfo) -> None:
        self.hMixedEvent = HistHandler.loadHist(mixedEventInfo)
        self.hMixedEvent.SetName('hMixed_invMass')

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
            n_bins_same = self.hSameEvent.GetNbinsX()
            sameEventIntegral = self.hSameEvent.Integral(1, n_bins_same+1)
            self.hSameEvent.Scale(1/sameEventIntegral)
        else:
            print(tc.GREEN+'[INFO]: '+tc.RESET+'No same event histogram provided')
        if self.hMixedEvent:
            n_bins_mixed = self.hMixedEvent.GetNbinsX()
            mixedEventIntegral = self.hMixedEvent.Integral(1, n_bins_mixed+1)
            self.hMixedEvent.Scale(1/mixedEventIntegral)
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
            sameEventIntegral = self.hSameEvent.Integral(low_bin, high_bin)
            mixedEventIntegral = self.hMixedEvent.Integral(low_bin, high_bin)
            self.hMixedEvent.Scale(sameEventIntegral/mixedEventIntegral)
        else:
            print(tc.GREEN+'[INFO]: '+tc.RESET+'No histogram provided')

    def rebin(self, rebin_factor:int=2) -> None:
        '''
            Rebin the sameEvent and the mixedEvent histograms.
        '''
        if self.hSameEvent:     self.hSameEvent.Rebin(rebin_factor)
        if self.hMixedEvent:    self.hMixedEvent.Rebin(rebin_factor)

    # TODO: implement cut variation methods
        
    def bkg_subtraction(self) -> None:
        '''
            Subtract the background from the signal.
        '''
        if not self.hSameEvent or not self.hMixedEvent:
            print(tc.RED+'[ERROR]: '+tc.RESET+'No histograms provided')
            return

        self.hSubtracted = self.hSameEvent.Clone()
        self.hSubtracted.SetName('hSubtracted_invMass')
        self.hSubtracted.Reset()

        for bin in range(1, self.hSameEvent.GetNbinsX()+1):
            sameValue = self.hSameEvent.GetBinContent(bin)
            mixedValue = self.hMixedEvent.GetBinContent(bin)
            sameError = self.hSameEvent.GetBinError(bin)
            mixedError = self.hMixedEvent.GetBinError(bin)
            self.hSubtracted.SetBinContent(bin, sameValue-mixedValue)
            self.hSubtracted.SetBinError(bin, np.sqrt(sameError**2 + mixedError**2))

    def save(self) -> None:
        '''
            Save the histograms in the output file.
        '''
        self.dir.cd()
        self.hSameEvent.Write()
        self.hMixedEvent.Write()
        self.hSubtracted.Write()

    def produce_plot(self, outFilePath:str) -> None:
        '''
            Produce a plot with the invariant mass distribution.
        '''
        
        raise NotImplementedError('Method not implemented yet')