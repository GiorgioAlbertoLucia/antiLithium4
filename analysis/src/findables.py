import numpy as np
from ROOT import TFile, TH1F, TH2F
import uproot

class Findables:

    def __init__(self, inFilePath:str, outFilePath:str):
        
        self.inFile = TFile.Open(inFilePath, 'read')
        self.outFile = TFile(outFilePath, 'recreate')
        print(f'Creating output file {outFilePath}...')

    def close(self):
        self.inFile.Close()
        self.outFile.Close()

    def __buildEfficiency(self, hPtTrue, hPtReco, name:str):

        if 'TH1' in str(type(hPtReco)):
            hEff = TH1F(name, f'{name}; p_{{T}} [GeV/#it{{c}}]; Efficiency', hPtTrue.GetNbinsX(), hPtTrue.GetXaxis().GetXmin(), hPtTrue.GetXaxis().GetXmax())
            for xbin in range(1, hPtTrue.GetNbinsX()):
                if hPtTrue.GetBinContent(xbin) > 0:
                    eff = hPtReco.GetBinContent(xbin)/hPtTrue.GetBinContent(xbin)
                    #effErr = np.sqrt(eff*(1-eff)/hPtTrue.GetBinContent(xbin))
                    hEff.SetBinContent(xbin, eff)
                    #hEff.SetBinError(xbin, effErr)
            return hEff

        if 'TH2' in str(type(hPtReco)):
            hEff = TH2F(name, f'{name}; p_{{T}} [GeV/#it{{c}}]; Efficiency', hPtTrue.GetNbinsX(), hPtTrue.GetXaxis().GetXmin(), hPtTrue.GetXaxis().GetXmax(), hPtReco.GetNbinsY(), hPtReco.GetYaxis().GetXmin(), hPtReco.GetYaxis().GetXmax())
            for ybin in range(1, hPtReco.GetNbinsY() + 1):   
                hEff.GetYaxis().SetBinLabel(ybin, hPtReco.GetYaxis().GetBinLabel(ybin)) 
                for xbin in range(1, hPtTrue.GetNbinsX()):
                    if hPtTrue.GetBinContent(xbin) > 0:
                        eff = hPtReco.GetBinContent(xbin, ybin)/hPtTrue.GetBinContent(xbin)
                        #effErr = np.sqrt(eff*(1-eff)/hPtTrue.GetBinContent(xbin))
                        hEff.SetBinContent(xbin, ybin, eff)
                        #hEff.SetBinError(xbin, ybin, effErr)
            return hEff

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

        inDir = self.inFile.Get('lithium4findables')

        for part in particles:
            hPtTrue = inDir.Get(f'ptTrue{part}')
            hPtReco = inDir.Get(f'ptReco{part}')
            hPtRecoSel = inDir.Get(f'ptRecoSel{part}')

            hPtTrues.append(hPtTrue)
            hPtRecos.append(hPtReco)
            hPtRecoSels.append(hPtRecoSel)
            hEfficiencies.append(self.__buildEfficiency(hPtTrue, hPtReco, f'efficiency{part}'))
            hEfficiencieSels.append(self.__buildEfficiency(hPtTrue, hPtRecoSel, f'efficiencySel{part}'))

        outDir.cd()
        for true, reco, eff, recoSel, effSel in zip(hPtTrues, hPtRecos, hEfficiencies, hPtRecoSels, hEfficiencieSels):
            true.Write()
            reco.Write()
            eff.Write() 
            recoSel.Write()
            effSel.Write()

    

