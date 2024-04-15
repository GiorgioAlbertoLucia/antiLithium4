from preprocessing import Preprocessor
from extraVisual import ExtraVisual

if __name__ == '__main__':

    inputFile = '/Users/glucia/Projects/ALICE/antiLithium4/MCWorkflowAnalysis/AO2D_lit_mc.root'
    cfgVisualFile = '/Users/glucia/Projects/ALICE/antiLithium4/MCanalysis/src/config/visualCfg.yml'
    preprocessor = Preprocessor()
    preprocessor.open(inputFile, treeName='O2lithium4tablemc', dirPrefix='DF*')
    preprocessor.defineVariables()
    preprocessor.visualize(cfgVisualFile)

    extraVisual = ExtraVisual(preprocessor.genDataset, preprocessor.recoDataset, cfgVisualFile)
    extraVisual.betheBloch()
    extraVisual.PIDinTracking()
    extraVisual.ptResolutionRec()
