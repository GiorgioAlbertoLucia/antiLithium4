'''
    Script to determin a Bethe Bloch-like parametrisation for the cluster size distribution
'''

import os
import sys
import yaml
import ctypes
import numpy as np
import polars as pl
from ROOT import (TFile, TDirectory, TH1F, TH2F, TF1, TCanvas, gInterpreter, 
                  TGraphErrors, TMultiGraph, kRed, kGreen, kBlue, kOrange, kCyan)
import logging
from typing import Dict, List, Tuple

# Include BetheBloch C++ file
CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
BETHEBLOCH_DIR = os.path.join(CURRENT_DIR, '..', 'include', 'BetheBloch.hh')
gInterpreter.ProcessLine(f'#include "{BETHEBLOCH_DIR}"')
from ROOT import BetheBloch, BetheBlochAleph

# Custom imports
sys.path.append('../..')
from framework.src.axis_spec import AxisSpec
from framework.src.hist_handler import HistHandler
from framework.src.graph_handler import GraphHandler
from framework.src.fitter import Fitter
from framework.utils.terminal_colors import TerminalColors as tc
from framework.utils.root_setter import obj_setter
from framework.utils.timeit import timeit
from utils.particles import ParticleMasses

lgging = logging.getLogger(__name__)

class BetheBlochParametrisation:

    def __init__(self, debug: bool = False):
        
        self.debug = debug
        self.config = None
        self.BetheBloch_params = {
            'kp1': -187.9255, 
            'kp2': -0.26878, 
            'kp3': 1.16252,
            'kp4': 1.15149, 
            'kp5': 2.47912
        }
        self.dir = None
        self.reset_fit_results()
        self.fitted_particle = None
        self.fitter = None

        self.h2 = None
        self.xs = None
        self.yerrs = None

    def load_config(self, config_file: str) -> None:
        if self.debug: print(tc.RED+'DEBUG:\t'+tc.RESET+'_load_config')
        with open(config_file, 'r') as file:
            self.config = yaml.safe_load(file)

    def init_config(self, cfg_label: str) -> None:
        if self.debug: print(tc.RED+'DEBUG:\t'+tc.RESET+'init_config')
        self.cfg = self.config[cfg_label]

    def _update_params(self, params: Dict[str, float]) -> None:
        if self.debug: print(tc.RED+'DEBUG:\t'+tc.RESET+'_update_params')
        for key, val in params.items():
            self.BetheBloch_params[key] = val

    def reset_fit_results(self) -> None:
        self.fit_results = pl.DataFrame({'x': pl.Series(values=[], dtype=pl.Float64), 
                                         'sx': pl.Series(values=[], dtype=pl.Float64), 
                                         'signal_norm': pl.Series(values=[], dtype=pl.Float64), 
                                         'signal_mean': pl.Series(values=[], dtype=pl.Float64), 
                                         'signal_sigma': pl.Series(values=[], dtype=pl.Float64),
                                         'signal_sigma_norm': pl.Series(values=[], dtype=pl.Float64),
                                         'fm_norm': pl.Series(values=[], dtype=pl.Float64),
                                         'fm_mean': pl.Series(values=[], dtype=pl.Float64), 
                                         'fm_sigma': pl.Series(values=[], dtype=pl.Float64),
                                         'fm_sigma_norm': pl.Series(values=[], dtype=pl.Float64),
                                         'fm_prob': pl.Series(values=[], dtype=pl.Float64),
                                         'fm_prob_err': pl.Series(values=[], dtype=pl.Float64),
                                         'chi2/ndf': pl.Series(values=[], dtype=pl.Float64)})

    def _create_axis_specs(self, cfg_plot: Dict[str, any], particle: str) -> Tuple[AxisSpec, AxisSpec]:
        ''' 
            Create axis specs for x and y variables from the configuration file
        '''
        axis_spec_x = AxisSpec(
            cfg_plot['nXBins'], cfg_plot['xMin'], cfg_plot['xMax'],
            f"{cfg_plot['name']}_{particle}", cfg_plot['title']
        )
        axis_spec_y = AxisSpec(
            cfg_plot['nYBins'], cfg_plot['yMin'], cfg_plot['yMax'],
            f"{cfg_plot['name']}_{particle}_y", f";{cfg_plot['yLabel']};counts"
        )
        return axis_spec_x, axis_spec_y
    
    def _set_output_dir(self, out_dir: TDirectory) -> None:
        if self.debug: print(tc.RED+'DEBUG:\t'+tc.RESET+'_set_output_dir')
        self.dir = out_dir

    def select_fit_particle(self, particle: str) -> None:
        if self.debug: print(tc.RED+'DEBUG:\t'+tc.RESET+'select_fit_particle')
        self.fitted_particle = particle

    ########### Fits

    def upload_h2(self, data) -> None:
        if self.debug: print(tc.RED+'DEBUG:\t'+tc.RESET+'upload_h2')

        if 'TH2' in str(type(data)):
            self.h2 = data
        elif 'DataFrame' in str(type(data)):
            hist_handler = HistHandler.createInstance(data)
            axis_spec_x, axis_spec_y = self._create_axis_specs(self.cfg['plot'], self.fitted_particle)
            self.h2 = hist_handler.buildTH2(self.cfg['plot']['xVariable'], self.cfg['plot']['yVariable'], axis_spec_x, axis_spec_y)

    def _get_fit_bins(self) -> Tuple[int, int, int]:
        if self.debug: print(tc.RED+'DEBUG:\t'+tc.RESET+'_get_fit_bins')
        first_bin = self.h2.GetXaxis().FindBin(self.cfg['xMinFit'])
        last_bin = self.h2.GetXaxis().FindBin(self.cfg['xMaxFit'])
        first_bin_double_fit = self.h2.GetXaxis().FindBin(self.cfg['xMinDoubleFit'])
        last_bin_double_fit = self.h2.GetXaxis().FindBin(self.cfg['xMaxDoubleFit'])
        return first_bin, last_bin, first_bin_double_fit, last_bin_double_fit

    def generate_bethe_bloch_points(self, **kwargs) -> None:
        '''
            Save graph to output file
        '''
        if self.debug: print(tc.RED+'DEBUG:\t'+tc.RESET+'generate_bethe_bloch_points')
        self._update_params(self.cfg['BBparams'])

        if self.h2 is None:
            raise ValueError('No histogram provided')

        first_bin, last_bin, first_bin_double_fit, last_bin_double_fit = self._get_fit_bins()
        for ibin in range(first_bin, last_bin + 1):
            h1 = self.h2.ProjectionY(f"slice_{ibin}", ibin, ibin)
            self._fit_slice(h1, ibin, **kwargs)
            del h1

        self.create_plots_from_fit_results()

    def _fit_slice(self, h1:TH1F, ibin: int, **kwargs) -> None:
        if self.debug: print(tc.RED+'DEBUG:\t'+tc.RESET+'_fit_histogram_slice')
        
        fit_status, fit = self._configure_fit_function(h1, ibin, **kwargs)
        
        canvas = TCanvas()
        h1.SetTitle(f"fit {self.h2.GetXaxis().GetBinCenter(ibin)}")
        h1.Draw("hist e0")

        if fit_status.IsValid():

            fit.SetLineColor(kRed)
            fit.Draw("same")
            if ibin < self.h2.GetXaxis().FindBin(self.cfg['xMaxDoubleFit']) and ibin >= self.h2.GetXaxis().FindBin(self.cfg['xMinDoubleFit']):
                sig, fm = self.create_signal_and_fake_matching_functions(fit)
                sig.SetLineColor(kOrange)
                sig.Draw("same")
                fm.SetLineColor(kGreen)
                fm.Draw("same")
                fm_prob, fm_prob_err = self.calculate_fake_matching_fraction(fit)

                fit_results = {
                    'x': self.h2.GetXaxis().GetBinCenter(ibin),
                    'sx': self.h2.GetXaxis().GetBinWidth(ibin)/2,
                    'signal_norm': fit.GetParameter(0),
                    'signal_mean': fit.GetParameter(1),
                    'signal_sigma': fit.GetParameter(2),
                    'signal_sigma_norm': fit.GetParameter(2) / np.sqrt(fit.GetParameter(0)),
                    'fm_norm': fit.GetParameter(self.cfg['signalFit']['nParams']),
                    'fm_mean': fit.GetParameter(self.cfg['signalFit']['nParams']+1),
                    'fm_sigma': fit.GetParameter(self.cfg['signalFit']['nParams']+2),
                    'fm_sigma_norm': fit.GetParameter(self.cfg['signalFit']['nParams']+2) / np.sqrt(fit.GetParameter(self.cfg['signalFit']['nParams'])),
                    'fm_prob': fm_prob,
                    'fm_prob_err': fm_prob_err,
                    'chi2/ndf': fit.GetChisquare() / fit.GetNDF()
                }
            else:
                fit_results = {
                    'x': self.h2.GetXaxis().GetBinCenter(ibin),
                    'sx': self.h2.GetXaxis().GetBinWidth(ibin)/2,
                    'signal_norm': fit.GetParameter(0),
                    'signal_mean': fit.GetParameter(1),
                    'signal_sigma': fit.GetParameter(2),
                    'signal_sigma_norm': fit.GetParameter(2) / np.sqrt(fit.GetParameter(0)),
                    'fm_norm': -1.,
                    'fm_mean': -1.,
                    'fm_sigma': -1.,
                    'fm_sigma_norm': -1.,
                    'fm_prob': -1.,
                    'fm_prob_err': -1.,
                    'chi2/ndf': fit.GetChisquare() / fit.GetNDF()
                }
            self.fit_results = pl.concat([self.fit_results, pl.DataFrame(fit_results)])
                
        self.dir.cd()
        canvas.SetTitle(f"fit_{self.h2.GetXaxis().GetBinCenter(ibin)}")
        canvas.BuildLegend()
        canvas.Write(f"fit_{self.h2.GetXaxis().GetBinCenter(ibin)}", )

    def create_plots_from_fit_results(self) -> None:
        if self.debug: print(tc.RED+'DEBUG:\t'+tc.RESET+'create_and_save_graphs')
        cfg_plot = self.cfg['plot']
        self.fit_results = self.fit_results.with_columns(signal_sigma_over_mean=(pl.col('signal_sigma') / pl.col('signal_mean')))
        graph_handler = GraphHandler(self.fit_results)

        signal_graph = graph_handler.createTGraphErrors('x', 'signal_mean', 'sx', 'signal_sigma_norm')
        obj_setter(signal_graph, name=f'{cfg_plot["name"]}_{self.fitted_particle}_sig_points', title=f'{cfg_plot["name"]}_{self.fitted_particle}_sig_points; {cfg_plot["xLabel"]}; {cfg_plot["yLabel"]}', marker_color=kRed, marker_style=20, marker_size=1)

        signal_sigma_over_mean_graph = graph_handler.createTGraph('x', 'signal_sigma_over_mean')
        obj_setter(signal_sigma_over_mean_graph, name=f'{cfg_plot["name"]}_{self.fitted_particle}_sig_sigma_points', title=f'{cfg_plot["name"]}_{self.fitted_particle}_sig_sigma_points; {cfg_plot["xLabel"]}; sigma / {cfg_plot["yLabel"]}', marker_color=kBlue, marker_style=20, marker_size=1)
        fit_sigma = TF1('fit_sigma', 'pol2')
        signal_sigma_over_mean_graph.Fit(fit_sigma, 'SM+')
        logging.info(f"Fit sigma: [0] = {fit_sigma.GetParameter(0)}")
        logging.info(f"Fit sigma: [1] = {fit_sigma.GetParameter(1)}")
        logging.info(f"Fit sigma: [2] = {fit_sigma.GetParameter(2)}")

        fit_results = self.fit_results.filter(pl.col('fm_mean') > 0)
        graph_handler = GraphHandler(fit_results)
        fm_graph = graph_handler.createTGraphErrors('x', 'fm_mean', 'sx', 'fm_sigma_norm')
        obj_setter(fm_graph, name=f'{cfg_plot["name"]}_{self.fitted_particle}_fm_points', title=f'{cfg_plot["name"]}_{self.fitted_particle}_fm_points; {cfg_plot["xLabel"]}; {cfg_plot["yLabel"]}', marker_color=kGreen, marker_style=20, marker_size=1)
        fm_prob_graph = graph_handler.createTGraphErrors('x', 'fm_prob', 'sx', 'fm_prob_err')
        obj_setter(fm_prob_graph, name=f'{cfg_plot["name"]}_{self.fitted_particle}_fm_prob_points', title=f'{cfg_plot["name"]}_{self.fitted_particle}_fm_prob_points; {cfg_plot["xLabel"]}; fake matching probability', marker_color=kOrange, marker_style=20, marker_size=1)

        self.dir.cd()
        signal_graph.Write(f"{cfg_plot['name']}_{self.fitted_particle}_sig_points")
        signal_sigma_over_mean_graph.Write(f"{cfg_plot['name']}_{self.fitted_particle}_sig_sigma_points")
        fm_graph.Write(f"{cfg_plot['name']}_{self.fitted_particle}_fm_points")
        fm_prob_graph.Write(f"{cfg_plot['name']}_{self.fitted_particle}_fm_prob_points")

    ########### Fits
        
    def _configure_fit_function(self, h1: TH1F, ibin: int, **kwargs) -> TF1:
        if self.debug: print(tc.RED+'DEBUG:\t'+tc.RESET+'_configure_fit_function')
        
        if ibin < self.h2.GetXaxis().FindBin(self.cfg['xMaxDoubleFit']) and ibin >= self.h2.GetXaxis().FindBin(self.cfg['xMinDoubleFit']):
            self.fitter = Fitter(h1, ['fmFit', 'signalFit'], self.cfg)
            self.fitter.auto_initialise()
        else:
            self.fitter = Fitter(h1, ['signalFit'], self.cfg)
            self.fitter.auto_initialise()
        
        if self.debug:  fit_status, fit = self.fitter.perform_fit(fit_option='RSLM+')
        else:           fit_status, fit = self.fitter.perform_fit(fit_option='RSLQM+')
        return fit_status, fit

    ########### Fake matching

    def create_signal_and_fake_matching_functions(self, fit: TF1) -> Tuple[TF1, TF1]:
        if self.debug: print(tc.RED+'DEBUG:\t'+tc.RESET+'create_signal_and_fake_matching_functions')
        sig = TF1('sig', 'gaus', self.cfg['yMinFit'], self.cfg['yMaxFit'])
        sig.SetParameters(fit.GetParameter(0), fit.GetParameter(1), fit.GetParameter(2))
        fm = TF1('fm', 'gaus', self.cfg['yMinFit'], self.cfg['yMaxFit'])
        offset_par = self.cfg['signalFit']['nParams']
        fm.SetParameters(fit.GetParameter(offset_par), fit.GetParameter(offset_par+1), fit.GetParameter(offset_par+2))
        return sig, fm

    def calculate_fake_matching_fraction(self, fit: TF1) -> Tuple[float, float]:
        '''
            Calculate the fake matching fraction as the ratio of the integral of the fake matching function

            OLD: Calculate the fake matching fraction as the ratio of the integral of the fake matching function 
            to the sum of the integrals of the signal and fake matching functions in the range of 2 sigmas
        '''
        if self.debug: print(tc.RED+'DEBUG:\t'+tc.RESET+'calculate_fake_matching_fraction')
        sig, fm = self.create_signal_and_fake_matching_functions(fit)
        #sig_center = sig.GetParameter(1)
        #sig_sigma = sig.GetParameter(2)
        #sig_frac = sig.Integral(sig_center - 2 * sig_sigma, sig_center + 2 * sig_sigma)
        #fm_frac = fm.Integral(sig_center - 2 * sig_sigma, sig_center + 2 * sig_sigma)

        sig_frac = sig.Integral(self.cfg['yMinFit'], self.cfg['yMaxFit'])
        fm_frac = fm.Integral(self.cfg['yMinFit'], self.cfg['yMaxFit'])
        fm_prob = fm_frac / (sig_frac + fm_frac)
        fm_prob_err = np.sqrt(fm_prob * (1 - fm_prob) / (sig_frac + fm_frac))
        if self.debug: print(tc.RED+'DEBUG:\t'+tc.RESET+f"Fake matching fraction: {fm_prob} +/- {fm_prob_err}")
        return fm_prob, fm_prob_err

    def save_fit_results(self, output_fit_path: str) -> None:
        if self.debug: print(tc.RED+'DEBUG:\t'+tc.RESET+'save_fit_results')
        self.fit_results.write_csv(output_fit_path)

    ########### Bethe-Bloch
          
    def _update_bethe_bloch_params(self, fit: TF1) -> None:
        if self.debug: print(tc.RED+'DEBUG:\t'+tc.RESET+'_update_bethe_bloch_params')
        logging.info("Updated Bethe-Bloch parameters:")
        for i, key in enumerate(self.BetheBloch_params.keys()):
            self.BetheBloch_params[key] = fit.GetParameter(key)
            logging.info(f"{key}: {fit.GetParameter(key)}")

    def upload_bethe_bloch_params(self, params: Dict[str, float]) -> None:
        if self.debug: print(tc.RED+'DEBUG:\t'+tc.RESET+'upload_bethe_bloch_params')
        self._update_params(params)

    def fit_bethe_bloch(self) -> None:
        '''
            Fit the TH2 by slices and fit the Bethe Bloch function to the mean values.
            Draw the Bethe Bloch fit and the distribution to a canvas.
        '''
        if self.debug: print(tc.RED+'DEBUG:\t'+tc.RESET+'fit_bethe_bloch')
        for key, val in self.cfg['BBparams'].items():
            self.BetheBloch_params[key] = val
            logging.info(f"{key}: {val}")
        
        bethe_bloch_fit = TF1('bethe_bloch', BetheBloch, self.cfg['xMinFit'], self.cfg['xMaxFit'], 5)
        bethe_bloch_fit.SetParameters(self.BetheBloch_params['kp1'], self.BetheBloch_params['kp2'], self.BetheBloch_params['kp3'], self.BetheBloch_params['kp4'], self.BetheBloch_params['kp5'])
        bethe_bloch_fit.SetParNames("kp1", "kp2", "kp3", "kp4", "kp5")

        mg = TMultiGraph()
        graph_handler = GraphHandler(self.fit_results)
        sig_points_graph = graph_handler.createTGraphErrors('x', 'signal_mean', 0, 'signal_sigma_norm')
        sig_points_graph.Fit(bethe_bloch_fit, 'RSM+')
        self._update_bethe_bloch_params(bethe_bloch_fit)

        self._plot_bethe_bloch_fit(sig_points_graph, mg, self.fitted_particle)
        self.draw_bethe_bloch_fit(self.fitted_particle)

        del bethe_bloch_fit, sig_points_graph, mg

    def _plot_bethe_bloch_fit(self, sig_points_graph: TGraphErrors, mg: TMultiGraph, particle: str) -> None:
        if self.debug: print(tc.RED+'DEBUG:\t'+tc.RESET+'_plot_bethe_bloch_fit')
        obj_setter(sig_points_graph, title=f"{particle} Bethe Bloch fit; {self.cfg['plot']['xLabel']}; {self.cfg['plot']['yLabel']}", marker_color=kRed, marker_style=20, marker_size=1)
        mg.Add(sig_points_graph)

        bethe_bloch_fit = TF1('bethe_bloch_fit', BetheBloch, self.cfg['xMinFit'], self.cfg['xMaxFit'], 5)
        bethe_bloch_fit.SetParameters(self.BetheBloch_params['kp1'], self.BetheBloch_params['kp2'], self.BetheBloch_params['kp3'], self.BetheBloch_params['kp4'], self.BetheBloch_params['kp5'])
        bethe_bloch_fit.SetParNames("kp1", "kp2", "kp3", "kp4", "kp5")

        canvas = TCanvas(f"cBB_{particle}", f"Bethe-Bloch fit for {particle}")
        mg.Draw("AP")
        bethe_bloch_fit.Draw("same")
        self.dir.cd()
        mg.Write(f"{particle}_BB_fit")

        del canvas, mg, bethe_bloch_fit

    def draw_bethe_bloch_fit(self, particle: str, **kwargs) -> None:
        if self.debug: print(tc.RED+'DEBUG:\t'+tc.RESET+'draw_bethe_bloch_fit')
        bethe_bloch_fit = TF1('bethe_bloch_fit', BetheBloch, self.cfg['xMinFit'], self.cfg['xMaxFit'], 5)
        bethe_bloch_fit.SetParameters(self.BetheBloch_params['kp1'], self.BetheBloch_params['kp2'], self.BetheBloch_params['kp3'], self.BetheBloch_params['kp4'], self.BetheBloch_params['kp5'])
        bethe_bloch_fit.SetParNames("kp1", "kp2", "kp3", "kp4", "kp5")

        canvas_name = kwargs.get('canvas_name', f"h2_cBB_{particle}")
        xmin_fit = kwargs.get('xmin_fit', self.h2.GetXaxis().GetXmin())
        xmax_fit = kwargs.get('xmax_fit', self.h2.GetXaxis().GetXmax())
        h2_canvas = TCanvas(canvas_name, f"Bethe-Bloch fit for {particle}")
        self.h2.Draw("colz")
        bethe_bloch_fit.SetRange(xmin_fit, xmax_fit)
        bethe_bloch_fit.Draw("same")
        h2_canvas.SetLogz()
        self.dir.cd()
        self.h2.Write(f'h2_BB_{particle}') 
        h2_canvas.Write()

        del h2_canvas, bethe_bloch_fit
    
    ########### NSigma distribution

    def _BBfunc(self, x: float) -> float:
        if self.debug: print(tc.RED+'DEBUG:\t'+tc.RESET+'_BBfunc')
        logging.info(f"BBfunc: {self.BetheBloch_params}")
        x = np.abs(x)
        beta = x / np.sqrt(1 + x**2)
        aa = beta**self.BetheBloch_params['kp4']
        bb = (1/x)**self.BetheBloch_params['kp5']
        bb = np.log(self.BetheBloch_params['kp3'] + bb)
        return (self.BetheBloch_params['kp2'] - aa - bb) * self.BetheBloch_params['kp1'] / aa

    def _get_sigma(self, x: float) -> float:
        if self.debug: print(tc.RED+'DEBUG:\t'+tc.RESET+'_get_sigma')
        idx = np.abs(self.xs - np.abs(x)).argmin()
        return self.yerrs[idx]

    def add_nsigma_column(self, dataset: pl.DataFrame, **kwargs) -> pl.DataFrame:
        '''
            Add an nsigma column to the dataset
        '''
        if self.debug: print(tc.RED+'DEBUG:\t'+tc.RESET+'create_nsigma_distribution')
        cfg_plot = self.cfg['plot']
        col_name_exp = f"fExp{cfg_plot['yVariable'][1:]}{self.fitted_particle}"
        logging.info(f"Bethe Bloch function: {self.BetheBloch_params}")

        self.xs = self.fit_results['x'].to_numpy()
        self.yerrs = self.fit_results['signal_sigma'].to_numpy()
        x_variable = cfg_plot['xVariable']
        if kwargs.get('convert_to_beta', False):
            x_variable = f"fBetaGamma{self.fitted_particle}"
            dataset = dataset.with_columns((pl.col('fPAbs')/ParticleMasses[self.fitted_particle]).alias(x_variable))
        dataset = dataset.with_columns(self._BBfunc(pl.col(x_variable)).alias(col_name_exp))
        
        sigma_params = kwargs.get('sigma_params', self.cfg['sigma_params']) # pol2 parametrisation for sigma/xvar
        sigma_x_var = kwargs.get('sigma_x_var', self.cfg['sigma_params']['xVar'])
        dataset = dataset.with_columns((pl.col(col_name_exp)*(sigma_params['kp0']+sigma_params['kp1']*pl.col(sigma_x_var)+sigma_params['kp2']*pl.col(sigma_x_var)*pl.col(sigma_x_var))).alias(f"fSigma{cfg_plot['yVariable'][1:]}"))
        dataset = dataset.with_columns(((pl.col(cfg_plot['yVariable']) - pl.col(col_name_exp)) / pl.col(f"fSigma{cfg_plot['yVariable'][1:]}")).alias(f"fNSigma{self.fitted_particle}"))

        return dataset

    def draw_sigma_distribution(self, dataset:pl.DataFrame, particle: str) -> None:
        '''
            Draw the sigma distribution.
            NOTE: The dataset should contain the sigma column and only entries for the particle to draw
        '''
        if self.debug: print(tc.RED+'DEBUG:\t'+tc.RESET+'draw_sigma_distribution')
        cfg_plot = self.cfg['plot_sigma']
        axis_spec_x, axis_spec_y = self._create_axis_specs(cfg_plot, particle)
        logging.info(f"Selected data\n{dataset.columns}\n{dataset[[cfg_plot['xVariable'], cfg_plot['yVariable']]].describe()}")
        hist_handler = HistHandler.createInstance(dataset)
        hist = hist_handler.buildTH2(cfg_plot['xVariable'], cfg_plot['yVariable'], axis_spec_x, axis_spec_y)
        self.dir.cd()
        hist.Write(f"{cfg_plot['name']}_{particle}_sigma")
        del hist, hist_handler

    def draw_nsigma_distribution(self, dataset:pl.DataFrame, particle: str) -> None:
        '''
            Draw the nsigma distribution.
            NOTE: The dataset should contain the nsigma column and only entries for the particle to draw
        '''
        if self.debug: print(tc.RED+'DEBUG:\t'+tc.RESET+'draw_nsigma_distribution')
        cfg_plot = self.cfg['plot_nsigma']
        axis_spec_x, axis_spec_y = self._create_axis_specs(cfg_plot, particle)
        sigma_col = f"fSigma{self.cfg['plot']['yVariable'][1:]}"
        logging.info(f"Selected data\n{dataset.columns}\n{dataset[[cfg_plot['xVariable'], cfg_plot['yVariable'], sigma_col, self.cfg['plot']['yVariable']]].describe()}")
        hist_handler = HistHandler.createInstance(dataset)
        hist = hist_handler.buildTH2(cfg_plot['xVariable'], cfg_plot['yVariable'], axis_spec_x, axis_spec_y)
        self.dir.cd()
        hist.Write(f"{cfg_plot['name']}_{particle}_nsigma")
        del hist, hist_handler

    def draw_expected_values(self, dataset:pl.DataFrame, particle:str) -> None:
        if self.debug: print(tc.RED+'DEBUG:\t'+tc.RESET+'draw_expected_values')
        cfg_plot = self.cfg['plot']
        col_name_exp = f"fExp{cfg_plot['yVariable'][1:]}{self.fitted_particle}"
        axis_spec_x, axis_spec_y = self._create_axis_specs(cfg_plot, 'Pr')
        hist_handler = HistHandler.createInstance(dataset)
        hist = hist_handler.buildTH2(cfg_plot['xVariable'], col_name_exp, axis_spec_x, axis_spec_y)
        self.dir.cd()
        obj_setter(hist, name=f"{cfg_plot['name']}_expected_{particle}", title=f"{cfg_plot['name']}_expected; {cfg_plot['xLabel']}; {cfg_plot['yLabel']}") 
        hist.Write(f"{cfg_plot['name']}_expected")
        del hist, hist_handler

    def _compute_purity_efficiency_in_slice(self, data: Dict[str, Dict[str, pl.DataFrame]], xlow: float, xup: float, nsigma_cut: float, other_particle: str) -> Tuple[float, float, float, float]:
        ''''
            Compute purity and efficiency in a slice of the Bethe-Bloch curve.

            Parameters:
            data: Dict[str, Dict[str, pl.DataFrame]] - dictionary of datasets -> [predicted_label][true_label]
            xlow: float - lower bound of the slice
            xup: float - upper bound of the slice
            nsigma_cut: float - cut on the nsigma distribution

            Returns:
            purity: float - purity of the slice 
            purity_err: float - error on the purity
            efficiency: float - efficiency of the slice
            eff_err: float - error on the efficiency
        '''
        if self.debug: print(tc.RED+'DEBUG:\t'+tc.RESET+'_compute_purity_efficiency_in_slice')
        cfg_plot = self.cfg['plot']

        true_pos = len(data[self.fitted_particle][self.fitted_particle].filter((np.abs(pl.col(cfg_plot['xVariable'])) > xlow) & (np.abs(pl.col(cfg_plot['xVariable'])) < xup)))
        true_neg = len(data[other_particle][other_particle].filter((np.abs(pl.col(cfg_plot['xVariable'])) > xlow) & (np.abs(pl.col(cfg_plot['xVariable'])) < xup)))
        false_pos = len(data[self.fitted_particle][other_particle].filter((np.abs(pl.col(cfg_plot['xVariable'])) > xlow) & (np.abs(pl.col(cfg_plot['xVariable'])) < xup)))
        false_neg = len(data[other_particle][self.fitted_particle].filter((np.abs(pl.col(cfg_plot['xVariable'])) > xlow) & (np.abs(pl.col(cfg_plot['xVariable'])) < xup)))

        #purity = 0 if (true_pos+false_pos) == 0 else true_pos/(true_pos+false_pos)
        # normalised purity
        purity = 0. if ((true_pos+false_pos) == 0 or (true_pos+false_neg == 0) or (false_pos+true_neg == 0)) else true_pos/(true_pos+false_neg)/(true_pos/(true_pos+false_neg)+false_pos/(false_pos+true_neg))
        efficiency = 0. if (true_pos+false_neg) == 0 else true_pos/(true_pos+false_neg)
        purity_err = 0. if (true_pos+false_pos) == 0 else np.sqrt(purity*(1-purity)/(true_pos+false_pos))
        eff_err = 0. if (true_pos+false_neg) == 0 else np.sqrt(efficiency*(1-efficiency)/(true_pos+false_neg))

        return purity, purity_err, efficiency, eff_err

    def purity_efficiency(self, dataset:pl.DataFrame, nsigma_cut:float, **kwargs) -> None:
        if self.debug: print(tc.RED+'DEBUG:\t'+tc.RESET+'purity_efficiency')
        cfg_plot = self.cfg['plot']

        other_particle = kwargs.get('other_particle', 'Pi')

        step = (cfg_plot['xMax'] - cfg_plot['xMin'])/cfg_plot['nXBins']

        dss = {
            other_particle: dataset.filter(pl.col('fPartID') == self.config['species'].index(other_particle)),
            self.fitted_particle: dataset.filter(pl.col('fPartID') == self.config['species'].index(self.fitted_particle))
        }

        dss_pred = {
            other_particle: {key:ds.filter(np.abs(pl.col(f'fNSigma{self.fitted_particle}')) > nsigma_cut) for key, ds in dss.items()},
            self.fitted_particle: {key:ds.filter(np.abs(pl.col(f'fNSigma{self.fitted_particle}')) < nsigma_cut) for key, ds in dss.items()}
        }

        variables_df = pl.DataFrame({'x': pl.Series(values=[], dtype=pl.Float64), 
                                     'purity': pl.Series(values=[], dtype=pl.Float64), 
                                     'purity_err': pl.Series(values=[], dtype=pl.Float64), 
                                     'efficiency': pl.Series(values=[], dtype=pl.Float64), 
                                     'eff_err': pl.Series(values=[], dtype=pl.Float64)
                                    })
                                     
        for ix in np.arange(cfg_plot['xMin'], cfg_plot['xMax'], step):
            p, p_err, e, e_err = self._compute_purity_efficiency_in_slice(dss_pred, ix, ix+step, nsigma_cut, other_particle)
            variables_value = {
                'x': ix+step/2,
                'purity': p,
                'purity_err': p_err,
                'efficiency': e,
                'eff_err': e_err
            }
            variables_df = pl.concat([variables_df, pl.DataFrame(variables_value)])

        eff_pur = TMultiGraph(f'eff_pur_{cfg_plot["xLabel"]}', f'Efficiency & Purity {self.fitted_particle}; {cfg_plot["xLabel"]}; Purity and Efficiency')
        graph_handler = GraphHandler(variables_df)
        purity = graph_handler.createTGraphErrors('x', 'purity', 0, 'purity_err')
        obj_setter(purity, name=f'purity_{cfg_plot["xLabel"]}', title=f'Purity {self.fitted_particle}; {cfg_plot["xLabel"]}; Purity', marker_color=kCyan-3, marker_style=20, marker_size=1)
        efficiency = graph_handler.createTGraphErrors('x', 'efficiency', 0, 'eff_err')
        obj_setter(efficiency, name=f'efficiency_{cfg_plot["xLabel"]}', title=f'Efficiency {self.fitted_particle}; {cfg_plot["xLabel"]}; Efficiency', marker_color=kOrange, marker_style=20, marker_size=1)

        eff_pur.Add(purity)
        eff_pur.Add(efficiency)

        canvas = TCanvas(f'c_eff_pur', f'Efficiency & Purity {self.fitted_particle}')
        canvas.cd()
        eff_pur.Draw('AP')
        canvas.BuildLegend()
        self.dir.cd()
        efficiency.Write(f'efficiency_{cfg_plot["xLabel"]}')
        purity.Write(f'purity_{cfg_plot["xLabel"]}')
        canvas.Write()



    ########### General
    
    @timeit
    def run_all(self, cfg_label: str, output_dir: TDirectory, fit_particle: str) -> None:
        if self.debug: print(tc.RED+'DEBUG:\t'+tc.RESET+'run_all')
        self._set_output_dir(output_dir)
        
        self.select_fit_particle(fit_particle)
        self.init_config(cfg_label)
        self.generate_bethe_bloch_points()
        self.fit_bethe_bloch()
        self.create_nsigma_distribution(self.fitted_particle)
        for part in self.config['species']:
            self.draw_nsigma_distribution(part)
        self.draw_expected_values()
        self.purity_efficiency(2)
            

