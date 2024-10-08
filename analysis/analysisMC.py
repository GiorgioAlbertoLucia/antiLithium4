import sys
sys.path.append('..')
from framework.src.data_handler import TableHandler, TaskHandler

from src.preprocessing import MCPreprocessor
from analysis.studies.studies import * 
from src.findables import Findables



def preprocessing(inFilePath, cfgVisualFile, correctH3:bool=True) -> MCPreprocessor:

    dataHandler = TableHandler(inFilePath=inFilePath, treeName='O2lithium4tablemc', dirPrefix='DF*')
    preprocessor = MCPreprocessor(dataHandler)
    preprocessor.defineVariables()
    preprocessor.defineKstar()
    if correctH3:   preprocessor.correctPtH3hp()
    preprocessor.visualize(cfgVisualFile)

    return preprocessor

def studies(preprocessor, cfgVisualFile) -> None:

    pidForTrkStudy = PIDforTrkStudy(preprocessor, cfgVisualFile)
    pidForTrkStudy.pidVsPtRes()

    ptResolutionStudy = PtResolutionStudy(preprocessor, cfgVisualFile)
    ptResolutionStudy.ptResolution()

    betheBlochStudy = BetheBlochStudy(preprocessor, cfgVisualFile)
    betheBlochStudy.fitBetheBloch()
    betheBlochStudy.drawBetheBloch()

    h3inTrkStudy = H3inTrkStudy(preprocessor, cfgVisualFile)
    h3inTrkStudy.fitH3pt()
    h3inTrkStudy.drawH3pt()

    Study.close()

def findablesStudies(inFilePath: str, outFilePath: str):

    dataHandler = TaskHandler(inFilePath, mainDir='lithium4findables')
    findables = Findables(dataHandler, outFilePath)
    findables.evaluateEfficiency()
    findables.close()


if __name__ == '__main__':

    print()
    #inFilePath = '/data/galucia/lithium4/MC/processed_AO2D_LHC24f3.root'
    inFilePath = '/Users/glucia/Projects/ALICE/antiLithium4/data/01/AO2D_processed.root'
    cfgVisualFile = '/Users/glucia/Projects/ALICE/antiLithium4/analysis/src/config/cfgMC.yml'

    preprocessor = preprocessing(inFilePath, cfgVisualFile, correctH3=False)
    
    studies(preprocessor, cfgVisualFile)

    #inputFileFindables = '/home/galucia/antiLithium4/task/MCWorkflowFindables/AnalysisResults_primaryP.root'
    #outputFileFindables = '/home/galucia/antiLithium4/analysis/output/mc_findables_primaryP.root'
    inputFileFindables = '/home/galucia/antiLithium4/task/MCWorkflowFindables/AnalysisResults_li4.root'
    outputFileFindables = '/home/galucia/antiLithium4/analysis/output/mc_findables_li4.root'
    #findablesStudies(inputFileFindables, outputFileFindables)
