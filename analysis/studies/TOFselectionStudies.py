'''
    Class to visualize different variables with particular selections on the dataset
'''
import numpy as np
from ROOT import TH2F, TF1, TObjArray

import sys
sys.path.append('..')
from .studies import StandaloneStudy

sys.path.append('../..')
from framework.src.hist_info import HistLoadInfo
from framework.src.hist_handler import HistHandler
from framework.utils.terminal_colors import TerminalColors as tc

class TOFselectionStudy(StandaloneStudy):

    def __init__(self, config, h2_mass_info:HistLoadInfo, particle:str='Pr'):

        super().__init__(config)
        self.dir = TOFselectionStudy.outFile_shared.mkdir('TOFselection')

        # Bethe-Bloch parameters
        
        # Parameters from data
        self.h2 = HistHandler.loadHist(h2_mass_info)
        self.particle = particle

        self.mean_points = None
        self.sigma_points = None
        self.resolution = None
        self.res_params = None
        
    def rebinx(self, rebin_factor:int=2) -> None:
        self.h2.RebinX(rebin_factor)

    def fit_slices(self, **kwargs) -> None:

        gaus = TF1('gaus', 'gaus')
        results = TObjArray()
        self.h2.FitSlicesY(gaus, 0, -1, 0, 'QL', results)

        self.mean_points = results[1]
        self.sigma_points = results[2]

        self.resolution = self.sigma_points.Clone('resolution'+self.particle)
        for ibin in range(1, self.resolution.GetNbinsX()+1):
            if self.mean_points.GetBinContent(ibin) == 0:
                self.resolution.SetBinContent(ibin, 0)
                continue
            self.resolution.SetBinContent(ibin, self.sigma_points.GetBinContent(ibin)/self.mean_points.GetBinContent(ibin))

        self.mean_points.SetTitle('; #it{p}_{T} (GeV/#it{c}); m_{TOF} (GeV/#it{c}^2)')
        self.sigma_points.SetTitle('; #it{p}_{T} (GeV/#it{c}); #sigma_{m_{TOF}} (GeV/#it{c}^2)')
        self.resolution.SetTitle('; #it{p}_{T} (GeV/#it{c}); #sigma_{m_{TOF}}/#mu_{m_{TOF}}')
    
    def fit_resolution(self) -> None:
        #self.res_fit = TF1('resolution_fit'+self.particle, 'pol2', self.resolution.GetXaxis().GetBinLowEdge(1), self.resolution.GetXaxis().GetBinLowEdge(self.resolution.GetNbinsX()+1))
        self.res_fit = TF1('resolution_fit'+self.particle, '[0]*exp([1]*abs(x))', self.resolution.GetXaxis().GetBinLowEdge(1), self.resolution.GetXaxis().GetBinLowEdge(self.resolution.GetNbinsX()+1))
        self.res_fit.SetParameter(0, 0.1)
        self.res_fit.SetParameter(1, 0.4)
        self.resolution.Fit(self.res_fit, 'RMS+')
        self.res_params = [self.res_fit.GetParameter(iparam) for iparam in range(2)]


    def fit_mass(self) -> None:
        if self.res_params:
            for ibin in range(1, self.mean_points.GetNbinsX()+1):
                counts = self.h2.ProjectionY('proj', ibin, ibin).GetEntries() if self.h2.ProjectionY('proj', ibin, ibin).GetEntries() > 0 else 1
                self.mean_points.SetBinError(ibin, self.mean_points.GetBinContent(ibin) * self.res_params[0]*np.exp(self.res_params[1]*abs(self.mean_points.GetBinCenter(ibin)) / counts))
        self.mass_fit = TF1('mass_fit'+self.particle, 'pol0', -3, 0)
        self.mean_points.Fit(self.mass_fit, 'RMS+') 

    def save(self) -> None:

        self.dir.cd()
        self.h2.Write()
        self.mean_points.Write('mean_TOFmass'+self.particle)
        self.sigma_points.Write('sigma_TOFmass'+self.particle)
        self.resolution.Write()
        self.mass_fit.Write()
        self.res_fit.Write()

