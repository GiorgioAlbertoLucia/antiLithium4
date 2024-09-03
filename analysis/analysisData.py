import sys
sys.path.append('..')
from framework.src.dataHandler import TableHandler, TaskHandler

from src.preprocessing import DataPreprocessor
from src.studies import * 
from src.invMassStudies import InvariantMassStudy
from src.findables import Findables



def preprocessing(inFilePath, cfgVisualFile, antimatterOnly=False) -> DataPreprocessor:

    dataHandler = TableHandler(inFilePath=inFilePath, treeName='O2lithium4table', dirPrefix='DF*')
    preprocessor = DataPreprocessor(dataHandler)
    preprocessor.defineVariables()
    if antimatterOnly: preprocessor.filterAntimatter()
    preprocessor.visualize(cfgVisualFile)

    return preprocessor

def studies(preprocessor, cfgVisualFile) -> None:

    betheBlochStudy = BetheBlochStudy(preprocessor, cfgVisualFile)
    betheBlochStudy.fitBetheBloch(rebin=2)
    betheBlochStudy.drawBetheBloch()

    invMassStudy = InvariantMassStudy(preprocessor, cfgVisualFile)
    invMassStudy.generalSelections()
    invMassStudy.invariantMass()    

    Study.close()

def findablesStudies(inFilePath: str, outFilePath: str):

    dataHandler = TaskHandler(inFilePath, mainDir='lithium4findables')
    findables = Findables(dataHandler, outFilePath)
    findables.evaluateEfficiency()
    findables.close()


if __name__ == '__main__':

    print()
    inFilePath = '/data/galucia/lithium4/same/LHC22o_pass4_minBias_Thin.root'
    cfgVisualFile = '/home/galucia/antiLithium4/analysis/src/config/cfgData.yml'

    preprocessor = preprocessing(inFilePath, cfgVisualFile, antimatterOnly=False)
    
    studies(preprocessor, cfgVisualFile)

    inputFileFindables = '/home/galucia/antiLithium4/task/MCWorkflowFindables/AnalysisResults.root'
    outputFileFindables = '/home/galucia/antiLithium4/analysis/output/MCfindables.root'
    #findablesStudies(inputFileFindables, outputFileFindables)
