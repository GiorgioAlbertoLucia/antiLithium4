import sys
sys.path.append('..')
from framework.src.dataHandler import TableHandler, TaskHandler

from src.preprocessing import DataPreprocessor
from src.studies import * 
from src.findables import Findables



def preprocessing(inFilePath, cfgVisualFile) -> DataPreprocessor:

    dataHandler = TableHandler(inFilePath=inFilePath, treeName='O2lithium4table', dirPrefix='DF*')
    preprocessor = DataPreprocessor(dataHandler)
    preprocessor.defineVariables()
    preprocessor.visualize(cfgVisualFile)

    return preprocessor

def studies(preprocessor, cfgVisualFile) -> None:

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
    inFilePath = '/data/galucia/lithium4/LHC23_pass4_skimmed_QC1/event_mixing/AO2D.root'
    cfgVisualFile = '/home/galucia/antiLithium4/analysis/src/config/cfgEventMixing.yml'

    preprocessor = preprocessing(inFilePath, cfgVisualFile)
    
    studies(preprocessor, cfgVisualFile)

    inputFileFindables = '/home/galucia/antiLithium4/task/MCWorkflowFindables/AnalysisResults.root'
    outputFileFindables = '/home/galucia/antiLithium4/analysis/output/MCfindables.root'
    #findablesStudies(inputFileFindables, outputFileFindables)
