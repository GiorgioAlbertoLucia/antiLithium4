'''
    Classes for correlation function studies
'''
import numpy as np
from ROOT import TH1F, TH2F, TCanvas, TLine, TLegend
from ROOT import kOrange, kGray

from .studies import StandaloneStudy

from torchic.core.histogram import load_hist
from torchic import HistLoadInfo
from torchic.utils.terminal_colors import TerminalColors as tc

import sys
sys.path.append('../..')
from utils.root_setter import obj_setter

class CorrelationStudy(StandaloneStudy):

    def __init__(self, config=None, outputFile=None, sameEvent=None, mixedEvent=None, opt: str = 'Anti'):
        '''
            Study to investigate the invariant mass distribution with different cuts.
        '''
        super().__init__(config, outputFile)
        if self.outFile:
            self.dir = self.outFile.mkdir(f'Correlation{opt}')

        self.h2SameEvent = None
        self.h2MixedEvent = None
        self.hSameEvent = None
        self.hMixedEvent = None

        if sameEvent:
            self.set_same_event(sameEvent)
        else:
            print(tc.MAGENTA+'[WARNING]: '+tc.RESET+'No same event histogram provided')
        if mixedEvent:
            self.set_mixed_event(mixedEvent)
        else:
            print(tc.MAGENTA+'[WARNING]: '+tc.RESET+'No mixedEvent histogram provided')

        self.hCorrelation = None
        self.hPull = None
        self.hGenuineCorrelation = None

        self.hSameEventCent = []
        self.hMixedEventCent = []
        self.hCorrelationCent = []
        self.hPullCent = []
        self.hGenuineCorrelationCent = []
        self.hCorrelationMassT = []

        self.CentralityBinEdges = [0, 10, 30, 50]
        self.MassTBinEdges = [3.747, 3.847, 3.947, 4.047]

    def clone_same_event(self, sameEvent) -> None:
        if 'TH2F' in str(type(sameEvent)):      self.h2SameEvent = sameEvent.Clone('h2Same_kstar')
        elif 'TH1F' in str(type(sameEvent)):    self.hSameEvent = sameEvent.Clone('hSame_kstar')

    def load_same_event(self, sameEventInfo:HistLoadInfo) -> None:
        print(tc.GREEN+'[INFO]: '+tc.RESET+f'Loading {sameEventInfo.hist_file_path}:{sameEventInfo.hist_name}')
        hist = load_hist(sameEventInfo)
        if 'TH2F' in str(type(hist)):                   self.h2SameEvent = hist.Clone('h2Same_kstar')
        elif 'TH1F' in str(type(hist)):                 self.hSameEvent = hist.Clone('hSame_kstar')

    def clone_mixed_event(self, mixedEvent) -> None:
        if 'TH2F' in str(type(mixedEvent)):     self.h2MixedEvent = mixedEvent.Clone('h2Mixed_kstar')
        elif 'TH1F' in str(type(mixedEvent)):   self.hMixedEvent = mixedEvent.Clone('hMixed_kstar')

    def load_mixed_event(self, mixedEventInfo:HistLoadInfo) -> None:
        print(tc.GREEN+'[INFO]: '+tc.RESET+f'Loading {mixedEventInfo.hist_file_path}:{mixedEventInfo.hist_name}')
        hist = load_hist(mixedEventInfo)
        if 'TH2F' in str(type(hist)):                   self.h2MixedEvent = hist.Clone('h2Mixed_kstar')
        elif 'TH1F' in str(type(hist)):                 self.hMixedEvent = hist.Clone('hMixed_kstar')

    def set_same_event(self, sameEvent) -> None:
        if ('TH2F' in str(type(sameEvent))) or ('TH1F' in str(type(sameEvent))):              
                                                        self.clone_same_event(sameEvent)
        elif 'HistLoadInfo' in str(type(sameEvent)):    self.load_same_event(sameEvent)
        else:                                           raise ValueError('Type not supported')
        if self.h2SameEvent:                            self.hSameEvent = self.h2SameEvent.ProjectionY('hSame_kstar')
    
    def set_mixed_event(self, mixedEvent) -> None:
        if ('TH2F' in str(type(mixedEvent))) or ('TH1F' in str(type(mixedEvent))):
                                                        self.clone_mixed_event(mixedEvent)
        elif 'HistLoadInfo' in str(type(mixedEvent)):   self.load_mixed_event(mixedEvent)
        else:                                           raise ValueError('Type not supported')
        if self.h2MixedEvent:                           self.hMixedEvent = self.h2MixedEvent.ProjectionY('hMixed_kstar')

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

    def normalize(self, low=0.5, high=0.9) -> float:
        '''
            Normalize the sameEvent and the mixedEvent histograms.
        '''

        if self.hMixedEvent and self.hSameEvent:
            return self._normalize_routine(self.hMixedEvent, self.hSameEvent, low, high)
        else:
            print(tc.GREEN+'[INFO]: '+tc.RESET+'No histogram provided')

    def _normalize_routine(self, histToNormalize:TH1F, histReference:TH1F, low:float, high:float) -> float:
        '''
            Normalize the sameEvent and the mixedEvent histograms.
        '''

        low_bin = histToNormalize.FindBin(low)
        low_edge = histToNormalize.GetBinLowEdge(low_bin)
        high_bin = histToNormalize.FindBin(high)
        high_edge = histToNormalize.GetBinLowEdge(high_bin+1)
        histIntegral = histToNormalize.Integral(low_bin, high_bin, 'width')
        referenceIntegral = histReference.Integral(low_bin, high_bin, 'width')
        if histIntegral < 1e-12:
            print(tc.RED+'[ERROR]: '+tc.RESET+'Normalization failed - denominator integral is zero')
            return 0.0
        if referenceIntegral < 1e-12:
            print(tc.RED+'[ERROR]: '+tc.RESET+'Normalization failed - numerator integral is zero')
            return 0.0

        histToNormalize.Scale(referenceIntegral/histIntegral)

        return referenceIntegral/histIntegral

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
            self.hSameEvent = self._custom_binning_routine(self.hSameEvent, bin_edges)
            
        if self.hMixedEvent:
            self.hMixedEvent = self._custom_binning_routine(self.hMixedEvent, bin_edges)

    def _custom_binning_routine(self, hist:TH1F, bin_edges:np.ndarray) -> TH1F:
        '''
            Define a custom binning for the histograms
        '''
        tmp_hist = TH1F(hist.GetName()+'_tmp', hist.GetTitle(), len(bin_edges)-1, bin_edges)
        for ibin in range(1, hist.GetNbinsX()+1):
            tmp_hist.Fill(hist.GetBinCenter(ibin), hist.GetBinContent(ibin))
        for ibin in range(1, tmp_hist.GetNbinsX()+1):
            tmp_hist.SetBinError(ibin, np.sqrt(tmp_hist.GetBinContent(ibin)))
        hist = tmp_hist.Clone(hist.GetName())
        del tmp_hist
        return hist

    def correlation_function(self) -> TH1F | None:
        '''
            Define the correlation function as the ratio bin by bin of the same event and the event mixing.
        '''
        if not self.hSameEvent or not self.hMixedEvent:
            print(tc.RED+'[ERROR]: '+tc.RESET+'No histograms provided')
            return
        
        self.hCorrelation = self._correlation_routine(self.hSameEvent, self.hMixedEvent)
        return self.hCorrelation

    def correlation_function_centrality(self, low_value_norm:float=0.5, high_value_norm:float=0.9, bin_edges:np.ndarray=[]) -> None:
        '''
            Define the correlation function as the ratio bin by bin of the same event and the event mixing.
        '''
        if not self.hSameEvent or not self.hMixedEvent:
            print(tc.RED+'[ERROR]: '+tc.RESET+'No histograms provided')
            return
        
        for i, cent_bin in enumerate(self.CentralityBinEdges[:-1]):
            sameEvent = self.h2SameEvent.ProjectionY(f'hSame_kstar_{cent_bin}', self.h2SameEvent.GetXaxis().FindBin(cent_bin), self.h2SameEvent.GetXaxis().FindBin(self.CentralityBinEdges[i+1]))
            mixedEvent = self.h2MixedEvent.ProjectionY(f'hMixed_kstar_{cent_bin}', self.h2MixedEvent.GetXaxis().FindBin(cent_bin), self.h2MixedEvent.GetXaxis().FindBin(self.CentralityBinEdges[i+1]))
            if len(bin_edges) > 0:
                sameEvent = self._custom_binning_routine(sameEvent, bin_edges)
                mixedEvent = self._custom_binning_routine(mixedEvent, bin_edges)
            self._normalize_routine(mixedEvent, sameEvent, low_value_norm, high_value_norm)
            self.hSameEventCent.append(sameEvent)
            self.hMixedEventCent.append(mixedEvent)
            correlation = self._correlation_routine(sameEvent, mixedEvent, suffix=f'_cent{self.CentralityBinEdges[i]}_{self.CentralityBinEdges[i+1]}')
            self.hCorrelationCent.append(correlation)

        self.hSameEvent = self.hSameEventCent[0].Clone('hSame_kstar')
        for i in range(1, len(self.hSameEventCent)):
            self.hSameEvent.Add(self.hSameEventCent[i])
        self.hMixedEvent = self.hMixedEventCent[0].Clone('hMixed_kstar')
        for i in range(1, len(self.hMixedEventCent)):
            self.hMixedEvent.Add(self.hMixedEventCent[i])
        
        self.hCorrelation = self._correlation_routine(self.hSameEvent, self.hMixedEvent, suffix='')

    def correlation_function_massT(self, sameEvent:HistLoadInfo, mixedEvent:HistLoadInfo, low_value_norm:float=0.5, high_value_norm:float=0.9, bin_edges:np.ndarray=[]) -> None:
        '''
            Define the correlation function as the ratio bin by bin of the same event and the event mixing.
        '''
        
        h2SameEvent = load_hist(sameEvent)
        h2MixedEvent = load_hist(mixedEvent)
        
        for i, massT_bin in enumerate(self.MassTBinEdges[:-1]):
            sameEvent = h2SameEvent.ProjectionY(f'hSame_kstar_{massT_bin}', h2SameEvent.GetXaxis().FindBin(massT_bin), h2SameEvent.GetXaxis().FindBin(self.MassTBinEdges[i+1]))
            mixedEvent = h2MixedEvent.ProjectionY(f'hMixed_kstar_{massT_bin}', h2MixedEvent.GetXaxis().FindBin(massT_bin), h2MixedEvent.GetXaxis().FindBin(self.MassTBinEdges[i+1]))
            if len(bin_edges) > 0:
                sameEvent = self._custom_binning_routine(sameEvent, bin_edges)
                mixedEvent = self._custom_binning_routine(mixedEvent, bin_edges)
            self._normalize_routine(mixedEvent, sameEvent, low_value_norm, high_value_norm)
            correlation = self._correlation_routine(sameEvent, mixedEvent, suffix=f'cent{self.MassTBinEdges[i]}_{self.MassTBinEdges[i+1]}')
            self.hCorrelationMassT.append(correlation)

    def _correlation_routine(self, sameEvent:TH1F, mixedEvent:TH1F, suffix:str='') -> TH1F:
        '''
            Define the correlation function as the ratio bin by bin of the same event and the event mixing.
        '''
        hCorrelation = sameEvent.Clone('hCorrelation_kstar'+suffix)
        hCorrelation.Reset()
        hCorrelation.SetTitle('k* Correlation; k* (GeV/#it{c}); C(k*)')

        for ibin in range(1, sameEvent.GetNbinsX()+1):
            valueSame = sameEvent.GetBinContent(ibin)
            valueMixed = mixedEvent.GetBinContent(ibin)
            errorSame = sameEvent.GetBinError(ibin)
            errorMixed = mixedEvent.GetBinError(ibin)
            
            if valueSame < 1e-12 or valueMixed < 1e-12: continue
            valueCorrelation = valueSame/valueMixed
            errorCorrelation = valueCorrelation*np.sqrt((errorSame/valueSame)*(errorSame/valueSame) + (errorMixed/valueMixed)*(errorMixed/valueMixed))
            hCorrelation.SetBinContent(ibin, valueCorrelation)
            hCorrelation.SetBinError(ibin, errorCorrelation)
        
        return hCorrelation

    def pull_distribution(self, genuine_correlation_load: HistLoadInfo = None) -> None:
        '''
            Compute the pull distribution
        '''
        if not self.hCorrelation:
            print(tc.RED+'[ERROR]: '+tc.RESET+'No correlation function defined')
            return

        self.hPull = self.hCorrelation.Clone('hPull_kstar')
        self.hPull.Reset()
        self.hPull.SetTitle('Pull distribution; k* (GeV/#it{c}); Pull')

        if genuine_correlation_load:
            self.hGenuineCorrelation = load_hist(genuine_correlation_load)
        else:
            self.hGenuineCorrelation = None

        for ibin in range(1, self.hCorrelation.GetNbinsX()+1):
            bin_center = self.hCorrelation.GetBinCenter(ibin)
            value = self.hCorrelation.GetBinContent(ibin)
            error = self.hCorrelation.GetBinError(ibin)
            if value < 1e-12: continue
            value_genuine = self.hGenuineCorrelation.GetBinContent(self.hGenuineCorrelation.FindBin(bin_center)) if self.hGenuineCorrelation else 1.0
            #error_genuine = self.hGenuineCorrelation.GetBinError(ibin) if self.hGenuineCorrelation else 0.0
            error_genuine = 0. if self.hGenuineCorrelation else 0.0
            self.hPull.SetBinContent(ibin, (value-value_genuine)/np.sqrt(error**2 + error_genuine**2))

    def pull_distribution_centrality(self, genuine_correlation_loads: list) -> None:

        for icent in range(len(self.CentralityBinEdges)-1):
            hPullCent = self.hCorrelationCent[icent].Clone(f'hPullCent_kstar_cent{self.CentralityBinEdges[icent]}_{self.CentralityBinEdges[icent+1]}')
            hPullCent.Reset()
            hPullCent.SetTitle(f'Pull distribution Centrality {self.CentralityBinEdges[icent]}-{self.CentralityBinEdges[icent+1]}; k* (GeV/#it{{c}}); Pull')

            hGenuineCorrelation = load_hist(genuine_correlation_loads[icent])
            hGenuineCorrelation.SetName(f'hGenuineCorrelationCent_kstar_cent{self.CentralityBinEdges[icent]}_{self.CentralityBinEdges[icent+1]}')
            self.hGenuineCorrelationCent.append(hGenuineCorrelation)

            for ibin in range(1, self.hCorrelationCent[icent].GetNbinsX()+1):
                bin_center = self.hCorrelationCent[icent].GetBinCenter(ibin)
                value = self.hCorrelationCent[icent].GetBinContent(ibin)
                error = self.hCorrelationCent[icent].GetBinError(ibin)
                if value < 1e-12: continue
                value_genuine = hGenuineCorrelation.GetBinContent(hGenuineCorrelation.FindBin(bin_center))
                #error_genuine = hGenuineCorrelation.GetBinError(ibin)
                error_genuine = 0.
                hPullCent.SetBinContent(ibin, (value-value_genuine)/np.sqrt(error**2 + error_genuine**2))
            self.hPullCent.append(hPullCent)


    def save(self, suffix='') -> None:

        self.dir.cd()
        if self.hSameEvent:             self.hSameEvent.Write(self.hSameEvent.GetName()+suffix)
        if self.hMixedEvent:            self.hMixedEvent.Write(self.hMixedEvent.GetName()+suffix)
        if self.hCorrelation:           self.hCorrelation.Write(self.hCorrelation.GetName()+suffix)
        if self.hPull:                  self.hPull.Write(self.hPull.GetName()+suffix)
        if self.hGenuineCorrelation:    self.hGenuineCorrelation.Write(self.hGenuineCorrelation.GetName()+suffix)
        for hCorrelationCent in self.hCorrelationCent:
            if hCorrelationCent:        hCorrelationCent.Write(hCorrelationCent.GetName()+suffix)
        for hPullCent in self.hPullCent:
            if hPullCent:               hPullCent.Write(hPullCent.GetName()+suffix)
        for hGenuineCorrelationCent in self.hGenuineCorrelationCent:
            if hGenuineCorrelationCent: hGenuineCorrelationCent.Write(hGenuineCorrelationCent.GetName()+suffix)
        for hCorrelationMassT in self.hCorrelationMassT:
            if hCorrelationMassT:       hCorrelationMassT.Write(hCorrelationMassT.GetName()+suffix)
        canvas = TCanvas('canvas', 'canvas', 800, 600)
        canvas.cd()
        
        
        if self.hMixedEvent:         
            self.hMixedEvent.SetMarkerColor(kGray+2)
            self.hMixedEvent.Draw('same hist')
        if self.hSameEvent:     
            self.hSameEvent.SetMarkerColor(kOrange-3)
            self.hSameEvent.Draw('same')

        self.dir.cd()
        canvas.Write('canvas'+suffix)

    def produce_plot(self, output_pdf:str) -> None:
        '''
            Produce the plots
        '''

        canvas = TCanvas('canvas', 'canvas', 800, 600)
        canvas.cd()
        hframe = canvas.DrawFrame(0.0, 0., 2.0, 7., ';k* (GeV/#it{c}); C(k*)')
        #hframe = canvas.DrawFrame(0.0, 0., 1.0, 2., ';k* (GeV/#it{c}); C(k*)')
        obj_setter(self.hCorrelation, marker_style=20, marker_color=797)
        self.hCorrelation.Draw('e1 same')
        const_line = TLine(0.0, 1.0, 2.0, 1.0)
        obj_setter(const_line, line_color=kGray+2, line_style=2, line_width=2)
        const_line.Draw('same')
        legend = TLegend(0.15, 0.6, 0.49, 0.8)
        legend.AddEntry(self.hCorrelation, 'p -^{3}He', 'p')
        legend.AddEntry(const_line, 'C(k*) = 1', 'l')
        legend.SetBorderSize(0)
        legend.SetTextSize(0.04)
        legend.Draw('same')
        canvas.SaveAs(output_pdf)
        print(tc.GREEN+'[INFO]: '+tc.RESET+'Correlation function plot saved in '+tc.UNDERLINE+tc.CYAN+f'{output_pdf}'+tc.RESET)