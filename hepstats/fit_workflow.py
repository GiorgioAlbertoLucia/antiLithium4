import os
import numpy as np
from copy import deepcopy
import zfit

import uproot
from hist import Hist

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import mplhep as hep

from scipy.stats import poisson

from core.utils import hist_to_tensor, pdf_to_numpy, space_from_hist, ratio_hist, pdf_to_hist

from hepstats.hypotests.calculators import AsymptoticCalculator
from hepstats.hypotests import UpperLimit, Discovery
from hepstats.hypotests.parameters import POI, POIarray

from torchic.utils.terminal_colors import TerminalColors as tc

class FitWorkflowCF:

    def __init__(self,  outfile: uproot.WritableDirectory, pdf_outpath: str):
        
        self.outfile = outfile
        self.outdir = outfile
        self.pdfpage = PdfPages(pdf_outpath)
        
        self.h_correlation_function = None
        self.CF_data = None
        self.CF_data_significance = None

        # zfit
        self.kstar = zfit.Space('kstar', limits=(0.02, 0.78))
        
        self.sig_counts = zfit.Parameter('sig_counts', 0.05, 0., 300)
        self.signal_pdf = None
        self.sig_params = {}
        self.sig_params_hesse = {} # Hesse errors
        self.bkg_counts = zfit.Parameter('bkg_counts', 1., 0., 1e6)
        self.bkg_datahist = None
        self.bkg_pdf = None
        self.model = None

        os.environ["ZFIT_DISABLE_TF_WARNINGS"] = "1"
        zfit.settings.advanced_warnings.all = False
        zfit.settings.changed_warnings.all = False

    def _log_decorator(message: str):
        def inner(func):
            def wrapper(*args, **kwargs):
                print(f'\n{tc.BG_GREEN}FitWorkflowCF:{tc.RESET} {message}')
                return func(*args, **kwargs)
            return wrapper
        return inner

    def add_directory(self, name: str):
        '''
            Create a directory in the output file
        '''
        self.outdir = name

    def load_hist_CF(self, h_correlation_function: Hist):
        '''
            Load the correlation function histogram
        '''
        
        self.h_correlation_function = h_correlation_function
        self.h_correlation_function.axes.name = ['kstar']
        self.kstar = space_from_hist(self.h_correlation_function, name='kstar')
        self.CF_data = zfit.data.BinnedData.from_hist(self.h_correlation_function)
        self.outfile[f'{self.outdir}/CF_data'] = self.CF_data.to_hist()

    @_log_decorator('Prefit the signal')
    def prefit_signal(self, h_signal: Hist, h_mixed_event: Hist):
        '''
            Fit the signal peak with a Crystal Ball function
        '''

        h_CFres = h_signal / h_mixed_event
        h_CFres.axes.name = ['prefit_kstar']   # it will be used as the space name for the data (needs to be consistent with the pdf)
        CFres_data = zfit.data.BinnedData.from_hist(h_CFres)
        space = space_from_hist(h_CFres, name='prefit_kstar')

        self.sig_params = {
            'sig_mean': zfit.Parameter('sig_mean', 0.0, -0.5, 0.5, label='#mu'),
            'sig_sigma': zfit.Parameter('sig_sigma', 0.1, 0., 1., label='#sigma'),
            'sig_aL': zfit.Parameter('sig_aL', 1.3, 0., 10., label='a_{L}'),
            'sig_nL': zfit.Parameter('sig_nL', 2.7, 0., 10., label='n_{L}'),
            'sig_aR': zfit.Parameter('sig_aR', 1.1, 0., 10., label='a_{R}'),
            'sig_nR': zfit.Parameter('sig_nR', 6., 0., 10., label='n_{R}'),
        }
        signal_pdf = self._build_signal_pdf(name='prefit_signal_pdf', obs=space)
        
        nll = zfit.loss.ExtendedBinnedNLL(model=signal_pdf.to_binned(space), data=CFres_data)
        minimizer = zfit.minimize.Minuit()
        result = minimizer.minimize(nll)
        self.sig_params_hesse = result.hesse()

        print(result)
        self.outfile[f'{self.outdir}/CFres_data'] = CFres_data.to_hist()
        hist = pdf_to_hist(signal_pdf.to_unbinned(), self.kstar.lower, self.kstar.upper)
        self.outfile[f'{self.outdir}/prefit_signal_pdf'] = hist
    
    @_log_decorator('Prepare the background template')
    def prepare_bkg_template(self, h_bkg: Hist):
        '''
            Get the Coulomb template for the p-He3 correlation function
        '''

        self._h_bkg_template = deepcopy(h_bkg)
        self._h_bkg_template.axes.name = ['kstar'] # it will be used as the space name for the pdf (needs to be consistent with the rest)
        self.bkg_pdf = zfit.pdf.HistogramPDF(self._h_bkg_template, 
                                             extended=self.bkg_counts,
                                             norm=None, name='bkg_pdf')
        self.outfile[f'{self.outdir}/bkg_pdf'] = self.bkg_pdf.to_hist()

    def pull_from_bkg(self, h_bkg):
        '''
            Pull the Coulomb template from the histogram
        '''

        return

    def _build_signal_pdf(self, name:str='signal_pdf', obs=None):

        if not obs:
            obs = self.kstar
        return zfit.pdf.DoubleCB(self.sig_params['sig_mean'], self.sig_params['sig_sigma'],
                                 self.sig_params['sig_aL'], self.sig_params['sig_nL'],
                                 self.sig_params['sig_aR'], self.sig_params['sig_nR'],
                                 name=name, obs=obs, 
                                 extended=self.sig_counts)
    
    def _init_signal_pdf(self):
        '''
            Initialize the signal pdf
        '''

        if not self.sig_params_hesse:
            raise ValueError('Signal parameters not set. Use prefit_signal() first.')
        
        mean = self.sig_params['sig_mean'].value()
        mean_err = self.sig_params_hesse[self.sig_params['sig_mean']]['error']
        mean_low, mean_up = mean - 3*mean_err, mean + 3*mean_err

        sigma = self.sig_params['sig_sigma'].value()
        sigma_low, sigma_up = sigma, 3*sigma
        
        self.sig_params = {
            'sig_mean': zfit.Parameter('sig_mean', mean, mean_low, mean_up, label='#mu'),
            'sig_sigma': zfit.Parameter('sig_sigma', sigma, sigma_low, sigma_up, label='#sigma'),
            'sig_aL': zfit.Parameter('sig_aL', self.sig_params['sig_aL'].value(), label='a_{L}', floating=False),
            'sig_nL': zfit.Parameter('sig_nL', self.sig_params['sig_nL'].value(), label='n_{L}', floating=False),
            'sig_aR': zfit.Parameter('sig_aR', self.sig_params['sig_aR'].value(), label='a_{R}', floating=False),
            'sig_nR': zfit.Parameter('sig_nR', self.sig_params['sig_nR'].value(), label='n_{R}', floating=False),
        }

        self.signal_pdf = self._build_signal_pdf(name='signal_pdf')

    @_log_decorator('\tPrefit the background')
    def _prefit_bkg(self):
        '''
            Prefit the p-He3 correlation function in a region away from the signal.
            This is used to fix the background parameters.
        '''

        if self.bkg_pdf is None:
            raise ValueError("Background PDF not initialized. Call prepare_bkg_template() first.")

        h_correlation_function_sliced = self.h_correlation_function[0.4j:0.8j]
        data_sliced = zfit.data.BinnedData.from_hist(h_correlation_function_sliced)

        h_bkg_sliced = self._h_bkg_template[0.4j:0.8j]
        bkg_pdf_sliced = zfit.pdf.HistogramPDF(h_bkg_sliced,
                                              extended=self.bkg_counts,
                                              norm=None,
                                              name='bkg_pdf_slice')

        nll = zfit.loss.ExtendedBinnedNLL(model=bkg_pdf_sliced, data=data_sliced)
        minimizer = zfit.minimize.Minuit()
        result = minimizer.minimize(nll)

        print(result)

        integral_bkg_range = self.bkg_pdf.integrate(limits=(0.4, 0.8))
        integral_bkg_total = self.bkg_pdf.integrate(limits=(0., 0.8))
        fitted_bkg_counts = result.params['bkg_counts']['value']
        self.bkg_counts.set_value(fitted_bkg_counts * integral_bkg_total / integral_bkg_range)
        self.bkg_counts.floating = False

    @_log_decorator('\tPlot the fit')
    def _plot_fit(self):
        
        fig, ax = plt.subplots(figsize=(8, 6))
        
        self.h_correlation_function.plot1d(ax=ax, label='Data', histtype='errorbar')

        self.bkg_pdf.to_binned(self.kstar).to_hist()\
            .plot1d(ax=ax, label='Background', histtype='step', linestyle='--')
        
        plt_signal_pdf = self._build_signal_pdf(name='plt_signal_pdf')
        xs, signal_ys = pdf_to_numpy(plt_signal_pdf, self.kstar, nbins=1000)
        ax.plot(xs, signal_ys * self.sig_counts / (self.sig_counts+self.bkg_counts),
                label='Signal', color='red', linestyle='--')
        
        self.model.to_binned(self.kstar).to_hist()\
            .plot1d(ax=ax, label='Model', histtype='step')

        ax.set_xlabel(r'\\textit{k}^* (GeV/\\textit{c})')
        ax.set_ylabel(r'C(\\textit{k}^*)')
        ax.legend()

        self.pdfpage.savefig(fig)

    def _plot_residuals(self):
        '''
            Plot the residuals of the fit
        '''

        return

    @_log_decorator('Fit the p-He3 correlation function')
    def fit_CF(self, get_upper_limit: bool = False, 
               get_pnull_significance: bool = False,
               get_pnull_significance_per_bin: bool = False):
        '''
            Fit the p-He3 correlation function
        '''

        if not self.CF_data:
            raise ValueError('Correlation function data not loaded. Use load_hist_CF() first.')
        
        self._init_signal_pdf()
        self._prefit_bkg()
        
        self.model = zfit.pdf.SumPDF(pdfs=[self.signal_pdf, self.bkg_pdf],
                                           #extended=self.sig_counts+self.bkg_counts,
                                           obs=self.kstar, name='model')

        nll = zfit.loss.ExtendedBinnedNLL(model=self.model, data=self.CF_data)
        minimizer = zfit.minimize.Minuit()
        result = minimizer.minimize(nll)
        print(result)

        self._plot_fit()
        if get_upper_limit:
            self._evaluate_upper_limit(nll, minimizer)
        if get_pnull_significance:
            self._evaluate_pnull_significance(nll, minimizer)
        if get_pnull_significance_per_bin:
            self._evaluate_local_significance()


    def integral_CF(self):
        '''
            Compute the integral of the k*^2 C(k*)
            This is used to compute the yield of Li4
        '''

        return

    @_log_decorator('\tEvaluate the upper limit')
    def _evaluate_upper_limit(self, loss, minimizer):

        calculator = AsymptoticCalculator(loss, minimizer, asimov_bins=100)
        poinull = POIarray(self.sig_counts, np.linspace(0.0, 25, 20))
        poialt = POI(self.sig_counts, 0)
        upper_limit = UpperLimit(calculator, poinull, poialt)
        print(f'{upper_limit.upperlimit(alpha=0.05, CLs=True)}')

    @_log_decorator('\tEvaluate the pnull significance')
    def _evaluate_pnull_significance(self, loss, minimizer):

        calculator = AsymptoticCalculator(loss, minimizer, asimov_bins=100)
        poinull = POI(self.sig_counts, 0.)
        discovery_test = Discovery(calculator, poinull)
        pnull, significance = discovery_test.result()
        print(f'{pnull=}, {significance=}')

    @_log_decorator('\tEvaluate the pnull significance per bin')
    def _evaluate_local_significance(self):
        '''
        Compute per-bin local significance comparing observed data to background-only prediction.
        Assumes bkg_counts has been fixed from prefit, and signal_pdf is active.
        '''

        data_hist = self.h_correlation_function
        data_vals = data_hist.values()

        original_sig_val = self.sig_counts.value()
        self.sig_counts.set_value(0.)
        self.sig_counts.floating = False

        model_bkg_only = zfit.pdf.SumPDF([self.signal_pdf, self.bkg_pdf], obs=self.kstar)
        hist_bkg_only = model_bkg_only.to_binned(self.kstar).to_hist()
        pred_vals = hist_bkg_only.values()

        self.sig_counts.set_value(original_sig_val)
        self.sig_counts.floating = True

        pvals = poisson.sf(data_vals - 1, mu=pred_vals)
        zvals = np.nan_to_num(poisson.isf(pvals, mu=pred_vals), nan=0.0, posinf=0.0)

        bin_centers = data_hist.axes[0].centers

        fig, ax = plt.subplots(figsize=(8, 6))
        ax.errorbar(bin_centers, zvals, fmt='o')
        ax.set_xlabel("k* (GeV/c)")
        ax.set_ylabel("Local significance [σ]")
        ax.set_title("Bin-by-bin signal excess significance")
        self.pdfpage.savefig(fig)
        
    

    def clear_data(self):
        '''
            Clear the data
        '''

        self.h_correlation_function = None
        self.CF_data = None

        self.signal_pdf = None
        self.sig_params = {}
        self.bkg_datahist = None
        self.bkg_pdf = None
        self.model = None

        self.pdfpage.close()

    def close_canvas(self):
        '''
            Close the canvas
        '''

        #self.pdf_canvas.Clear()
        #self.pdf_canvas.Print(f'{self.pdf_outpath})')
