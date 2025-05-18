import numpy as np
from hist import Hist

class EHist:
    def __init__(self, hist:Hist):
        self._hist = hist
        self._errors = None

    def set_errors(self, errors):
        if np.shape(errors) != self.hist.shape:
            raise ValueError(f"Shape mismatch: expected {self.hist.shape}, got {np.shape(errors)}")
        self._errors = np.array(errors)

    @property
    def hist(self):
        return self._hist

    @property
    def errors(self):
        return self._errors

    def fill_from_root(self, hist):
        '''
            Fill the EHist from a ROOT histogram. Only supports 1D histograms.
        '''
        
        yerrs = []
        for ibin in range(1, hist.GetNbinsX() + 1):
            self._hist.fill(hist.GetBinCenter(ibin), weight=hist.GetBinContent(ibin))
            yerrs.append(hist.GetBinError(ibin))
        self.set_errors(yerrs)

            