from analysisFramework.src.dataHandler import TableHandler

from src.preprocessing import Preprocessor
from src.extraVisual import ExtraVisual
from src.findables import Findables

if __name__ == '__main__':

    print()
    inFilePath = '/home/galucia/antiLithium4/MCWorkflowAnalysis/AO2D_lit_mc.root'
    cfgVisualFile = '/home/galucia/antiLithium4/MCanalysis/src/config/visualCfg.yml'

    dataHandler = TableHandler(inFilePath=inFilePath, treeName='O2lithium4tablemc', dirPrefix='DF*')
    preprocessor = Preprocessor(dataHandler)
    preprocessor.defineVariables()
    preprocessor.visualize(cfgVisualFile)
#
    #extraVisual = ExtraVisual(preprocessor.genDataset, preprocessor.recoDataset, cfgVisualFile)
    #extraVisual.betheBloch()
    #extraVisual.PIDinTracking()
    #extraVisual.ptResolutionRec()
    #extraVisual.close()

    #inputFileFindables = '/home/galucia/antiLithium4/MCWorkflowFindables/AnalysisResults.root'
    #outputFileFindables = '/home/galucia/antiLithium4/MCanalysis/output/MCfindables.root'
    #findables = Findables(inputFileFindables, outputFileFindables)
    #findables.evaluateEfficiency()
    #findables.close()
