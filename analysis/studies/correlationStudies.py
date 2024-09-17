'''
    Classes for invariant mass studies
'''
from ROOT import TFile, TH1F

from .studies import Study

import sys
sys.path.append('..')
from ..src.preprocessing import Preprocessor

sys.path.append('../..')
from framework.src.axis_spec import AxisSpec
from framework.src.hist_handler import HistHandler
from framework.utils.terminal_colors import TerminalColors as tc

class CorrelationStudy(Study):

    def __init__(self, preprocessor: Preprocessor, config):
        '''
            Study to investigate the invariant mass distribution with different cuts.
        '''
        super().__init__(preprocessor, config)
        self.dir = CorrelationStudy.outFile_shared.mkdir('correlation')

        cfg = self.config['Kstar']
        self.axisSpecX = AxisSpec(cfg['nXBins'], cfg['xMin'], cfg['xMax'], cfg['name'], cfg['title'])
        self.hMixing = TH1F('KstarMixing', '; K* (GeV/#it{c}); Counts', self.axisSpecX.nbins, self.axisSpecX.xmin, self.axisSpecX.xmax)
        self.hSameEvent = TH1F('KstarSame', '; K* (GeV/#it{c}); Counts', self.axisSpecX.nbins, self.axisSpecX.xmin, self.axisSpecX.xmax)

    
    def normalizeEventMixingBkg(self, sameEventPath:str, sameEventName:str)  -> None:
        '''
            Normalize the event mixing background to the sameEvent.

            Args:
                sameEventPath: path to the sameEvent file
                sameEventName: name of the sameEvent histogram
                lowInvMass: lower limit of the invariant mass
                upperInvMass: upper limit of the invariant mass
        '''
        sameEventFile = TFile(sameEventPath, 'read')
        hSameEvent = sameEventFile.Get(sameEventName)
        sameEventIntegral = hSameEvent.Integral()

        for x in self.dataset['full']['fKstar']: self.hMixing.Fill(x)
        mixingIntegral = self.hMixing.Integral()

        self.hMixing.Scale(sameEventIntegral/mixingIntegral)

        self.dir.cd()
        self.hMixing.Write('KstarMixingNormalized')

    def correlationFunction(self, sameEventPath:str, sameEventName:str) -> None:
        '''
            Define the correlation function as the ratio bin by bin of the same event and the event mixing.

            Args:
                sameEventPath: path to the signal file
                sameEventName: name of the signal histogram
        '''
        sameEventFile = TFile(sameEventPath, 'read')
        hSameEvent = sameEventFile.Get(sameEventName)
        hMixing = self.dir.Get('KstarMixingNormalized')

        if hSameEvent.GetNbinsX() != hMixing.GetNbinsX():
            raise ValueError('Histograms have different number of bins')
        
        hCorrelation = hSameEvent.Clone('KstarCorrelation')
        hCorrelation.Divide(hMixing)

        self.dir.cd()
        hSameEvent.Write('KstarSame')
        hCorrelation.Write()
        