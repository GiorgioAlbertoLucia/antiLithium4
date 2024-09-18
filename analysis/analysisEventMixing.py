import sys
sys.path.append('..')
from framework.src.data_handler import TableHandler, TaskHandler

from src.preprocessing import DataPreprocessor


def preprocessing(inFilePath, cfgVisualFile, antimatterOnly=False) -> DataPreprocessor:

    dataHandler = TableHandler(inFilePath=inFilePath, treeName='O2lithium4table', dirPrefix='DF*')
    preprocessor = DataPreprocessor(dataHandler)
    preprocessor.defineVariables()
    preprocessor.defineKstar()
    preprocessor.visualize(cfgVisualFile)

    if antimatterOnly: preprocessor.filterAntimatter()
    preprocessor.selectionsHe3()
    preprocessor.visualize(cfgVisualFile, output_suffix='_selectionsHe3')

    return preprocessor

if __name__ == '__main__':

    print()
    #inFilePath = '/data/galucia/lithium4/EM/LHC22o_pass4_minBias_Thin.root'
    inFilePath = ['/Users/glucia/Projects/ALICE/data/lithium/mixing/LHC24af_pass1_skimmed_mixing.root',
                  '/Users/glucia/Projects/ALICE/data/lithium/mixing/LHC24aj_pass1_skimmed_mixing.root',
                  '/Users/glucia/Projects/ALICE/data/lithium/mixing/LHC24al_pass1_skimmed_mixing.root']
    #cfgVisualFile = '/home/galucia/antiLithium4/analysis/src/config/cfgEventMixing.yml'
    cfgVisualFile = '/Users/glucia/Projects/ALICE/antiLithium4/analysis/src/config/cfgEventMixing.yml'

    preprocessor = preprocessing(inFilePath, cfgVisualFile, antimatterOnly=True)

    inputFileFindables = '/home/galucia/antiLithium4/task/MCWorkflowFindables/AnalysisResults.root'
    outputFileFindables = '/home/galucia/antiLithium4/analysis/output/MCfindables.root'
    #findablesStudies(inputFileFindables, outputFileFindables)
