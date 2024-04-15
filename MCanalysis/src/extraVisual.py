'''
    Class to visualize different variables with particular selections on the dataset
'''

from ROOT import TFile, TH1F, TH2F, TF1, TCanvas, gInterpreter, TObjArray

gInterpreter.ProcessLine('#include "../include/BetheBloch.hh"')
from ROOT import BetheBloch

import os
import sys
sys.path.append('..')
from utils.particle import PID

# Bethe-Bloch parameters
'''
# Parameters from data
BetheBlochParams = {'kp1':  -170.929,
                    'kp2':  -0.02554,
                    'kp3':  1.58698,
                    'kp4':  0.97078,
                    'kp5':  2.91625,
                    #'resolution':  0.09
        }
'''
# Parameters from MC
BetheBlochParams = {'kp1':  -190.936,
                    'kp2':  -0.22594,
                    'kp3':  1.17667,
                    'kp4':  1.16131,
                    'kp5':  2.43627,
                    #'resolution':  0.09
        }

class ExtraVisual:

    def __init__(self, genData:pd.DataFrame, recoData:pd.DataFrame, config) -> None:
        '''
            - genData (pd.DataFrame):  generated data to be visualized

        '''
        
        self.genData = genData
        self.recoData = recoData
        with open(config, 'r') as file:     self.config = yaml.safe_load(file)

        outFilePath = os.path.splitext(self.config['outputFilePath'])[0]+'Extra'+os.path.splitext(self.config['outputFilePath'])[1]
        self.outFile = TFile(outFilePath, 'recreate')
        print(f'Creating output file {outFilePath}...')
        

        # private methods
    
    # private methods
        
    def __buildHist(self, variable, name, title, nBins, xMin, xMax, data=None) -> TH1F:
        '''
            Build and return a histogram with given specifics. If no data is provided, generated data is used. 
        '''
        
        if data is None:    data = self.genData

        hist = TH1F(name, title, nBins, xMin, xMax)
        for value in data[variable]:    hist.Fill(value)
        return hist
    
    def __buildHist2(self, xVariable, yVariable, name, title, nXBins, xMin, xMax, nYBins, yMin, yMax, data=None) -> TH2F:
        '''
            Build and return a histogram with given specifics. If no data is provided, generated data is used. 
        '''

        if data is None:    data = self.genData
        hist2 = TH2F(name, title, nXBins, xMin, xMax, nYBins, yMin, yMax)
        for x, y in zip(data[xVariable], data[yVariable]):    hist2.Fill(x, y)
        return hist2
    


    def PIDinTracking(self) -> None:

        dir = self.outFile.mkdir('PIDinTracking')
        dir.cd()

        cfg = self.config['ptRes']
        for part in ['He3', 'Pr']:
            for pidHP in np.delete(self.data[f'fPIDtrk{part}'].unique(), np.where(self.data[f'fPIDtrk{part}'].unique() == 0xFFFFF)):
                title = cfg['title'].split(';', 1)[0] + f' {part} - PID trk hp {PID[pidHP]["label"]};' + cfg['title'].split(';', 1)[1]
                data = self.data.query(f'fPIDtrk{part} == {pidHP}', inplace=False)
                hist2 = self.__buildHist2(cfg['xVariable']+f'{part}', cfg['yVariable']+f'{part}', cfg['name']+f'{part}_pid_{PID[pidHP]["label"]}', title, cfg['nXBins'], cfg['xMin'], cfg['xMax'], cfg['nYBins'], cfg['yMin'], cfg['yMax'], data)
                hist2.Write()

    def ptResolutionRec(self) -> None:
        '''
            pt resolution plots. Non reconstructed particles are ignored
        '''

        dir = self.outFile.mkdir('ptResolutionRec')
        dir.cd()
        
        cfg = self.config['ptRes']
        for part in cfg['particle']:
            title = cfg['title'].split(';', 1)[0] + f' {part} (Reconstructed only);' + cfg['title'].split(';', 1)[1]
            hist2 = self.__buildHist2(cfg['xVariable']+part, cfg['yVariable']+part, cfg['name']+part, title, cfg['nXBins'], cfg['xMin'], cfg['xMax'], cfg['nYBins'], cfg['yMin'], cfg['yMax'], self.recoData)
            hist2.Write()

    def betheBloch(self, option='fit') -> None:

        dir = self.outFile.mkdir('BetheBloch')
        cfg = self.config['BetheBloch']

        dEdx = self.__buildHist2(cfg['xVariable']+'He3', cfg['yVariable']+'He3', cfg['name']+f' He3', cfg['title'], cfg['nXBins'], cfg['xMin'], cfg['xMax'], cfg['nYBins'], cfg['yMin'], cfg['yMax'])
        BBcurve = TF1(f'BetheBlochHe3', BetheBloch, cfg['xMin'], cfg['xMax'], len(BetheBlochParams.values()))
        for i, (parName, param) in enumerate(BetheBlochParams.items()):    
            BBcurve.SetParameter(i, param)
            BBcurve.SetParName(i, parName)
        
        if option == 'fit':
            gaus = TF1('gaus', 'gaus', cfg['yMin'], cfg['yMax'])
            results = TObjArray()
            dEdx.FitSlicesY(gaus, 0, -1, 0, 'Q', results)
            
            BBhist = results[1]
            BBres = results[2]
            for i in range(1, BBhist.GetNbinsX()+1):    BBhist.SetBinError(i, BBres.GetBinContent(i))
            BBhist.Fit(BBcurve, 'RM+')

        canvas = TCanvas(f'BBHe3', f'Bethe Bloch curve - He3')
        hframe = canvas.DrawFrame(cfg['xMin'], cfg['yMin'], cfg['xMax'], cfg['yMax'], f'Bethe Bloch He3; #beta #gamma; #frac{{dE}}{{dX}} (a.u.)')
        for i, (parName, param) in enumerate(BetheBlochParams.items()):    print(f'{parName}: {BBcurve.GetParameter(i)}')

        dir.cd()    
        BBcurve.Write()
        dEdx.Write()
            
        canvas.cd()
        dEdx.Draw()
        BBcurve.Draw('same')
        dir.cd()
        canvas.Write()
        

        


