'''
    Classes for invariant mass studies
'''
import numpy as np

from ROOT import TH1F, TH2F, TCanvas, TLine, TBox, TLegend, TF1
from ROOT import kGray, kOrange, kRed, gStyle

from .studies import StandaloneStudy

from torchic import HistLoadInfo
from torchic.core.histogram import load_hist
from torchic.utils.terminal_colors import TerminalColors as tc

import sys
sys.path.append('../..')
from utils.root_setter import obj_setter

class InvariantMassStudy(StandaloneStudy):

    def __init__(self, config, outputFile, sameEvent=None, mixedEvent=None, **kwargs):
        '''
            Study to investigate the invariant mass distribution with different cuts.
        '''
        super().__init__(config, outputFile)
        self.opt = kwargs.get('opt', '')
        self.dir = self.outFile.mkdir('InvariantMass'+self.opt)

        self.h2SameEvent = None
        self.h2MixedEvent = None
        self.hSameEvent = None
        self.hMixedEvent = None

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

        self.hSameEventCent = []
        self.hMixedEventCent = []
        self.hCorrectionCent = []
        self.hSubtractedCent = []
        self.hPullCent = []
        self.hRatioCent = []

        self.CentralityBinEdges = [0, 10, 30, 50]

    def clone_same_event(self, sameEvent:TH2F) -> None:
        self.hSameEvent = sameEvent.Clone('h2Same_invMass')

    def load_same_event(self, sameEventInfo:HistLoadInfo) -> None:
        print(tc.GREEN+'[INFO]: '+tc.RESET+f'Loading {sameEventInfo.hist_file_path}:{sameEventInfo.hist_name}')
        hist = load_hist(sameEventInfo)
        if 'TH2' in str(type(hist)):                    self.h2SameEvent = hist.Clone('h2Same_invMass')
        elif 'TH1' in str(type(hist)):                  self.hSameEvent = hist.Clone('hSame_invMass')
        else:                                           raise ValueError('Type not supported')

    def clone_mixed_event(self, mixedEvent:TH2F) -> None:
        self.h2MixedEvent = mixedEvent.Clone('hMixed_invMass')

    def load_mixed_event(self, mixedEventInfo:HistLoadInfo) -> None:
        print(tc.GREEN+'[INFO]: '+tc.RESET+f'Loading {mixedEventInfo.hist_file_path}:{mixedEventInfo.hist_name}')
        hist = load_hist(mixedEventInfo)
        if 'TH2' in str(type(hist)):                    self.h2MixedEvent = hist.Clone('h2Mixed_invMass')
        elif 'TH1' in str(type(hist)):                  self.hMixedEvent = hist.Clone('hMixed_invMass')

    def set_same_event(self, sameEvent) -> None:
        if 'TH2' in str(type(sameEvent)):               self.clone_same_event(sameEvent)
        elif 'HistLoadInfo' in str(type(sameEvent)):    self.load_same_event(sameEvent)
        else:                                           raise ValueError('Type not supported')
        if self.h2SameEvent:                            self.hSameEvent = self.h2SameEvent.ProjectionY('hSame_invMass')
    
    def set_mixed_event(self, mixedEvent) -> None:
        if 'TH2' in str(type(mixedEvent)):              self.clone_mixed_event(mixedEvent)
        elif 'HistLoadInfo' in str(type(mixedEvent)):   self.load_mixed_event(mixedEvent)
        else:                                           raise ValueError('Type not supported')
        if self.h2MixedEvent:                           self.hMixedEvent = self.h2MixedEvent.ProjectionY('hMixed_invMass')

    def set_centrality_binning(self, bin_edges:np.ndarray) -> None:
        self.CentralityBinEdges = bin_edges

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

    def _normalize_routine(self, histToNormalize:TH1F, histReference:TH1F, low:float, high:float) -> float:
        '''
            Normalize a histogram to a reference histogram in a given range.
        '''
        low_bin = histReference.FindBin(low)
        low_edge = histReference.GetBinLowEdge(low_bin)
        high_bin = histReference.FindBin(high)
        high_edge = histReference.GetBinLowEdge(high_bin+1)
        histIntegral = histToNormalize.Integral(low_bin, high_bin, 'width')
        referenceIntegral = histReference.Integral(low_bin, high_bin, 'width')
        if histIntegral < 1e-12: 
            print(tc.RED+'[ERROR]: '+tc.RESET+'Normalization failed - denominator integral is zero')
            return 0.
        if referenceIntegral < 1e-12:
            print(tc.RED+'[ERROR]: '+tc.RESET+'Normalization failed - numerator integral is zero')
            return 0.
        
        histToNormalize.Scale(referenceIntegral/histIntegral)
        return referenceIntegral/histIntegral

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

    def set_corrections_centrality(self, hCorrections: list) -> None:
        '''
            Set the correction histograms for different centrality bins.
        '''
        for ihist, hist in enumerate(hCorrections):
            tmp_hist = hist.Clone(f'hCorrection_cent{self.CentralityBinEdges[ihist]}_{self.CentralityBinEdges[ihist+1]}')
            self.hCorrectionCent.append(tmp_hist)

    def _correct_mixed_event_routine(self, histMixedEvent:TH1F, histCorrection:TH1F) -> None:
        '''
            Correct the mixed event histogram with a correction histogram.
            The correction histogram takes into account the effect of the interaction (obtained from the correlation function).
        '''

        for ibin in range(1, histMixedEvent.GetNbinsX()+1):
            mixedValue = histMixedEvent.GetBinContent(ibin)
            correctionValue = histCorrection.GetBinContent(ibin)
            histMixedEvent.SetBinContent(ibin, mixedValue*correctionValue)
            histMixedEvent.SetBinError(ibin, correctionValue*np.sqrt(mixedValue))

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

    def _custom_binning_routine(self, hist:TH1F, bin_edges:np.ndarray) -> TH1F:
        '''
            Define a custom binning for the histograms
        '''
        tmp_hist = TH1F(hist.GetName()+'_tmp', hist.GetTitle()+'_tmp', len(bin_edges)-1, bin_edges)
        for ibin in range(1, hist.GetNbinsX()):
            tmp_hist.Fill(hist.GetBinCenter(ibin), hist.GetBinContent(ibin))
        for ibin in range(1, tmp_hist.GetNbinsX()+1):
            tmp_hist.SetBinError(ibin, np.sqrt(tmp_hist.GetBinContent(ibin)))
        return tmp_hist

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

    def _bkg_subtraction_routine(self, histSameEvent:TH1F, histMixedEvent:TH1F, suffix:str='') -> TH1F:
        '''
            Subtract the background from the signal.
        '''

        histSubtracted = histSameEvent.Clone()
        histSubtracted.SetName('hSubtracted'+self.opt+'_invMass'+suffix)
        histSubtracted.Reset()

        for bin in range(1, histSameEvent.GetNbinsX()+1):
            sameValue = histSameEvent.GetBinContent(bin)
            mixedValue = histMixedEvent.GetBinContent(bin)
            sameError = histSameEvent.GetBinError(bin)
            mixedError = histMixedEvent.GetBinError(bin)
            histSubtracted.SetBinContent(bin, sameValue-mixedValue)
            histSubtracted.SetBinError(bin, np.sqrt(sameError**2 + mixedError**2))

        return histSubtracted

    def pull_distribution(self, low_edge_fit:float=None, high_edge_fit:float=None) -> None:
        '''
            Calculate the Pull distribution.
        '''
        if not self.hSubtracted:
            print(tc.RED+'[ERROR]: '+tc.RESET+'No subtracted histogram provided')
            return

        if not low_edge_fit: low_edge_fit = self.hSubtracted.GetBinLowEdge(1)
        if not high_edge_fit: high_edge_fit = self.hSubtracted.GetBinLowEdge(self.hSubtracted.GetNbinsX())

        pol0 = TF1('pol0', 'pol0', low_edge_fit, high_edge_fit)
        pol0.SetParameter(0, 0)
        self.hSubtracted.Fit(pol0, 'RMN+')

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

    def _pull_distribution_routine(self, histSubtracted:TH1F, low_edge_fit:float=None, high_edge_fit:float=None, suffix:str='') -> TH1F:
        '''
            Calculate the Pull distribution.
        '''
        if not low_edge_fit: low_edge_fit = histSubtracted.GetBinLowEdge(1)
        if not high_edge_fit: high_edge_fit = histSubtracted.GetBinLowEdge(histSubtracted.GetNbinsX())

        pol0 = TF1('pol0', 'pol0', low_edge_fit, high_edge_fit)
        pol0.SetParameter(0, 0)
        histSubtracted.Fit(pol0, 'RMN+')

        param = pol0.GetParameter(0)
        error = pol0.GetParError(0)

        histPull = histSubtracted.Clone()
        histPull.SetName('hPull'+self.opt+'_invMass'+suffix)
        histPull.Reset()

        for ibin in range(1, histSubtracted.GetNbinsX()+1):
            subtract_value = histSubtracted.GetBinContent(ibin)
            subtract_error = histSubtracted.GetBinError(ibin)
            pull = (subtract_value-param)/np.sqrt(subtract_error**2 + error**2)
            histPull.SetBinContent(ibin, pull)
            histPull.SetBinError(ibin, 0)

        return histPull, param

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

    def _ratio_distribution_routine(self, histSameEvent:TH1F, histMixedEvent:TH1F, suffix:str='') -> TH1F:
        '''
            Calculate the ratio distribution.
        '''
        histRatio = histSameEvent.Clone()
        histRatio.SetName('hRatio'+self.opt+'_invMass'+suffix)
        histRatio.Reset()

        for ibin in range(1, histSameEvent.GetNbinsX()+1):
            sameValue = histSameEvent.GetBinContent(ibin)
            mixedValue = histMixedEvent.GetBinContent(ibin)
            ratio = sameValue/mixedValue if mixedValue != 0 else 0
            ratioError = ratio*np.sqrt((histSameEvent.GetBinError(ibin)/sameValue)**2 + (histMixedEvent.GetBinError(ibin)/mixedValue)**2) if mixedValue != 0 and sameValue != 0 else 0
            histRatio.SetBinContent(ibin, ratio)
            histRatio.SetBinError(ibin, ratioError)

        return histRatio

    def invariant_mass_centrality(self, low_value_norm:float=3.78, high_value_norm:float=3.84, bin_edges:np.ndarray=[], rebin_factor:int=None) -> None:
        '''
            Perform the invariant mass analysis for different centrality bins.
        '''
        
        for i, centBin in enumerate(self.CentralityBinEdges[:-1]):
            sameEvent = self.h2SameEvent.ProjectionY(f'hSame_invMass_{centBin}', self.h2SameEvent.GetXaxis().FindBin(centBin), self.h2SameEvent.GetXaxis().FindBin(self.CentralityBinEdges[i+1]))
            mixedEvent = self.h2MixedEvent.ProjectionY(f'hMixed_invMass_{centBin}', self.h2MixedEvent.GetXaxis().FindBin(centBin), self.h2MixedEvent.GetXaxis().FindBin(self.CentralityBinEdges[i+1]))
            print(tc.YELLOW+'[DEBUG]: '+tc.RESET+f'Centrality bin {centBin}-{self.CentralityBinEdges[i+1]}')
            print(tc.YELLOW+'[DEBUG]: '+tc.RESET+f'Same Event bins {self.h2SameEvent.GetXaxis().FindBin(centBin)}-{self.h2SameEvent.GetXaxis().FindBin(self.CentralityBinEdges[i+1])}')
            print(tc.YELLOW+'[DEBUG]: '+tc.RESET+f'Mixed Event bins {self.h2MixedEvent.GetXaxis().FindBin(centBin)}-{self.h2MixedEvent.GetXaxis().FindBin(self.CentralityBinEdges[i+1])}')
            
            if rebin_factor:
                sameEvent.Rebin(rebin_factor)
                mixedEvent.Rebin(rebin_factor)

            if len(self.hCorrectionCent) > 0:
                self._correct_mixed_event_routine(mixedEvent, self.hCorrectionCent[i])

            self._normalize_routine(mixedEvent, sameEvent, low_value_norm, high_value_norm)

            if len(bin_edges) > 0:
                sameEvent = self._custom_binning_routine(sameEvent, bin_edges)
                mixedEvent = self._custom_binning_routine(mixedEvent, bin_edges)
            

            subtraction = self._bkg_subtraction_routine(sameEvent, mixedEvent, suffix=f'_cent{self.CentralityBinEdges[i]}_{self.CentralityBinEdges[i+1]}')
            pull, pol0_param = self._pull_distribution_routine(subtraction, high_edge_fit=high_value_norm, suffix=f'_cent{self.CentralityBinEdges[i]}_{self.CentralityBinEdges[i+1]}')
            print(tc.GREEN+'[INFO]: '+tc.RESET+f'Centrality bin {centBin}-{self.CentralityBinEdges[i+1]} pol0 fit: {pol0_param}')
            ratio = self._ratio_distribution_routine(sameEvent, mixedEvent, suffix=f'_cent{self.CentralityBinEdges[i]}_{self.CentralityBinEdges[i+1]}')

            self.hSameEventCent.append(sameEvent)
            self.hMixedEventCent.append(mixedEvent)
            self.hSubtractedCent.append(subtraction)
            self.hPullCent.append(pull)
            self.hRatioCent.append(ratio)

    def save(self, suffix='') -> None:
        '''
            Save the histograms in the output file.
        '''
        self.dir.cd()

        if self.hSameEvent:     self.hSameEvent.Write(self.hSameEvent.GetName()+suffix)
        if self.hMixedEvent:     self.hMixedEvent.Write(self.hMixedEvent.GetName()+suffix)
        if self.hCorrection:    self.hCorrection.Write(self.hCorrection.GetName()+suffix)
        if self.hSubtracted:    self.hSubtracted.Write(self.hSubtracted.GetName()+suffix)
        if self.hPull:          self.hPull.Write(self.hPull.GetName()+suffix)
        if self.hRatio:         self.hRatio.Write(self.hRatio.GetName()+suffix)

        for hSameEventCent in self.hSameEventCent:      hSameEventCent.Write()
        for hMixedEventCent in self.hMixedEventCent:    hMixedEventCent.Write()

        for hCorrectionCent in self.hCorrectionCent:    hCorrectionCent.Write()
        for hSubtractedCent in self.hSubtractedCent:    hSubtractedCent.Write()
        for hPullCent in self.hPullCent:                hPullCent.Write()
        for hRatioCent in self.hRatioCent:              hRatioCent.Write()
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