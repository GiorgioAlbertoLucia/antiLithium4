import sys
sys.path.append('..')
from framework.src.data_handler import TableHandler, TaskHandler

from src.preprocessing import DataPreprocessor
from analysis.studies.studies import Study
from analysis.studies.betheBlochStudies import BetheBlochStudy
from analysis.studies.clusterStudies import ClusterSizeParamStudy
from analysis.studies.invMassStudies import InvariantMassStudy
from src.findables import Findables



def preprocessing(inFilePath, cfgVisualFile, antimatterOnly=False) -> DataPreprocessor:

    dataHandler = TableHandler(inFilePath=inFilePath, treeName='O2lithium4table', dirPrefix='DF*')
    preprocessor = DataPreprocessor(dataHandler)
    preprocessor.defineVariables()
    preprocessor.defineKstar()
    if antimatterOnly: preprocessor.filterAntimatter()
    preprocessor.visualize(cfgVisualFile)

    preprocessor.selectionsHe3()
    preprocessor.visualize(cfgVisualFile, output_suffix='_selectionsHe3')
    #preprocessor.selectionsPr()
    #preprocessor.visualize(cfgVisualFile, output_suffix='_selectionsPr')

    return preprocessor

def studies(preprocessor, cfgVisualFile, bbConfig) -> None:

    clusterStudy = ClusterSizeParamStudy(preprocessor, cfgVisualFile, bbConfig)
    clusterStudy.drawBetheBloch('OLD_PARAMETERS_h2_cBB_Pr')
    clusterStudy.fitBetheBloch()

    #betheBlochStudy = BetheBlochStudy(preprocessor, cfgVisualFile)
    #betheBlochStudy.fitBetheBloch(rebin=2)
    #betheBlochStudy.drawBetheBloch()
#
    #invMassStudy = InvariantMassStudy(preprocessor, cfgVisualFile)
    #invMassStudy.generalSelections()
    #invMassStudy.invariantMass()    

    Study.close()

def findablesStudies(inFilePath: str, outFilePath: str):

    dataHandler = TaskHandler(inFilePath, mainDir='lithium4findables')
    findables = Findables(dataHandler, outFilePath)
    findables.evaluateEfficiency()
    findables.close()


if __name__ == '__main__':

    print()
    #inFilePath = '/Users/glucia/Projects/ALICE/antiLithium4/MCWorkflowAnalysis/AO2D_lit_mc.root'
    #inFilePath = '/data/galucia/lithium4/same/LHC22o_pass4_minBias_Thin.root'
    inFilePath = ['/Users/glucia/Projects/ALICE/data/lithium/same/LHC24af_pass1_skimmed_same.root',
                  '/Users/glucia/Projects/ALICE/data/lithium/same/LHC24aj_pass1_skimmed_same.root',
                  '/Users/glucia/Projects/ALICE/data/lithium/same/LHC24al_pass1_skimmed_same.root']
    cfgVisualFile = '/Users/glucia/Projects/ALICE/antiLithium4/analysis/src/config/cfgData.yml'
    #cfgVisualFile = '/home/galucia/antiLithium4/analysis/src/config/cfgData.yml'
    cfgBBFile = '/Users/glucia/Projects/ALICE/antiLithium4/analysis/src/config/cfgBetheBloch.yml'

    preprocessor = preprocessing(inFilePath, cfgVisualFile, antimatterOnly=True)
    
    studies(preprocessor, cfgVisualFile, cfgBBFile)

    #inputFileFindables = '/home/galucia/antiLithium4/task/MCWorkflowFindables/AnalysisResults.root'
    #outputFileFindables = '/home/galucia/antiLithium4/analysis/output/MCfindables.root'
    #findablesStudies(inputFileFindables, outputFileFindables)
