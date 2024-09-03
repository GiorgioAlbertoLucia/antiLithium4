import numpy as np
from ROOT import TFile, TH1F, TH2F
import uproot

import sys
sys.path.append('../..')
from framework.src.dataHandler import DataHandler
from framework.src.histHandler import HistHandler
from framework.utils.terminalColors import TerminalColors as tc

class Findables:

    def __init__(self, dataHandler: DataHandler, outFilePath:str):
        
        self.dataHandler = dataHandler
        self.data = self.dataHandler.inData
        self.histHandler = HistHandler.createInstance(self.data)

        self.outFile = TFile(outFilePath, 'recreate')
        print('Creating output file '+tc.UNDERLINE+tc.CYAN+f'{outFilePath}'+tc.RESET+'...')

    def close(self):
        # self.inFile.Close()
        self.outFile.Close()

    
    def evaluateEfficiency(self):
        '''
            Comment from PWG-LF meeting: expected efficiency of He3 is ~60% (TPC only, with TOF it goes down to 35%)
        '''

        outDir = self.outFile.mkdir('efficiency')
        particles = ['He3', 'Pr']
        hPtTrues = []
        hPtRecos = []
        hPtRecoSels = []
        hEfficiencies = []
        hEfficiencieSels = []
        hEfficiencieSelNormalized = []

        for part in particles:
            hPtTrue = self.histHandler.buildTH1(f'ptTrue{part}')
            hPtReco = self.histHandler.buildTH1(f'ptReco{part}')
            hPtRecoSel = self.histHandler.buildTH1(f'ptRecoSel{part}')

            hPtTrues.append(hPtTrue)
            hPtRecos.append(hPtReco)
            hPtRecoSels.append(hPtRecoSel)
            hEfficiency = self.histHandler.buildEfficiency(hPtReco, hPtTrue)
            hEfficiency.SetTitle(f'efficiency{part};p_{{T}}^{part} (GeV/#it{{c}});Efficiency')
            hEfficiencies.append(hEfficiency)

            effSel = self.histHandler.buildEfficiency(hPtRecoSel, hPtTrue)
            labels = {idx-1: hPtRecoSel.GetYaxis().GetBinLabel(idx) for idx in range(1, hPtRecoSel.GetNbinsY()+1)}
            self.histHandler.setLabels(effSel, labels, 'y')
            effSel.SetTitle(f'efficiencySel{part};p_{{T}}^{part} (GeV/#it{{c}});Selections')
            hEfficiencieSels.append(effSel)

            effSelNormalized = effSel.Clone(f'efficiencySelNormalized{part}')
            effSelNormalized.SetTitle(f'efficiencySelNormalized{part};p_{{T}}^{part} (GeV/#it{{c}});Selections')
            for ix in range(1, effSelNormalized.GetNbinsX()+1):
                for iy in range(effSelNormalized.GetNbinsY(), 0, -1):
                    if effSel.GetBinContent(ix, 1) == 0: continue
                    effSelNormalized.SetBinContent(ix, iy, effSel.GetBinContent(ix, iy)/effSel.GetBinContent(ix, 1))
            hEfficiencieSelNormalized.append(effSelNormalized)

        outDir.cd()
        for true, reco, eff, recoSel, effSel, effSelNormalized in zip(hPtTrues, hPtRecos, hEfficiencies, hPtRecoSels, hEfficiencieSels, hEfficiencieSelNormalized):
            true.Write()
            reco.Write()
            eff.Write() 
            recoSel.Write()
            effSel.Write()
            effSelNormalized.Write()

    

