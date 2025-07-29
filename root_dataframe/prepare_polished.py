import argparse
import yaml
import numpy as np
import ROOT
from ROOT import TFile, TChain, TList, TTree

from utils.particles import ParticleMasses
from utils.histogram_registry import HistogramRegistry, RegistryEntry
from utils.histogram_archive import register_qa_histograms, register_kstar_histograms, register_kstar_matter_histograms, \
    register_kstar_antimatter_histograms

ROOT.gROOT.LoadMacro('Common.h++')
from ROOT import NsigmaTpcHe, NsigmaITSHe, NsigmaITSPr, NsigmaTOFPr, averageClusterSize, CorrectPidTrkHe, \
  Kstar, ExpectedClusterSizeCosLambdaHe, ExpectedClusterSizeCosLambdaPr

ROOT.EnableImplicitMT(30)
ROOT.gROOT.SetBatch(True)

#parser = argparse.ArgumentParser(description='Configure the parameters of the script.')
#parser.add_argument('--config-file', dest='config_file', help='path to the YAML file with configuration.', default='')
#args = parser.parse_args()
#if args.config_file == '':
#    print('** No config file provided. Exiting. **')
#    exit()

#with open(args.config_file, 'r') as stream:
#  confFile = yaml.safe_load(stream)
#  conf = confFile["data_preparation"]

config_file = 'config/config_prepare.yml'
config = yaml.safe_load(open(config_file, 'r'))

base_selection = '(fSignedPtHe3 > 0 && fSignedPtHad > 0) || (fSignedPtHe3 < 0 && fSignedPtHad < 0)'
selections = config['selections']

selection = selections[0]
for sel in selections[1:]:
    selection += (' && ' + sel)
print(f'Selection: {selection}')

############################################################################################################

input_data = config['input_data']
tree_name = config['tree_name']
mode = config['mode']

file_data_list = input_data if isinstance(input_data, list) else [input_data]
chainData = TChain("O2he3hadtable")

for fileName in file_data_list:
  fileData = TFile(fileName)

  if mode == 'DF':
    for key in fileData.GetListOfKeys():
      keyName = key.GetName()
      if 'DF_' in keyName :
          print(f'Adding {fileName}/{keyName}/{tree_name} to the chain')
          chainData.Add(f'{fileName}/{keyName}/{tree_name}')
  elif mode == 'tree':
    print(f'Adding {fileName}/{tree_name} to the chain')
    chainData.Add(f'{fileName}/{tree_name}')

############################################################################################################################################################################
rdf = ROOT.ROOT.RDataFrame(chainData) \
      .Define('fSignedPtHad', 'fPtHad') \
      .Define('fSignHe3', 'fPtHe3/std::abs(fPtHe3)') \
      .Redefine('fPtHe3', 'std::abs(fPtHe3)') \
      .Redefine('fPtHad', 'std::abs(fPtHad)') \
      .Redefine('fPtHe3', '(fPIDtrkHe3 == 7) || (fPtHe3 > 2.5) ? fPtHe3 : CorrectPidTrkHe(fPtHe3)') \
      .Define('fSignedPtHe3', 'fPtHe3 * fSignHe3') \
      .Define(f'fEHe3', f'std::sqrt((fPtHe3 * std::cosh(fEtaHe3))*(fPtHe3 * std::cosh(fEtaHe3)) + {ParticleMasses["He"]}*{ParticleMasses["He"]})') \
      .Define(f'fEHad', f'std::sqrt((fPtHad * std::cosh(fEtaHad))*(fPtHad * std::cosh(fEtaHad)) + {ParticleMasses["Pr"]}*{ParticleMasses["Pr"]})') \
      .Define('fDeltaEta', 'fEtaHe3 - fEtaHad') \
      .Define('fDeltaPhi', 'fPhiHe3 - fPhiHad') \
      .Redefine('fInnerParamTPCHe3', 'fInnerParamTPCHe3 * 2') \
      .Define('fNSigmaTOFHad', 'NsigmaTOFPr(fPtHad * std::cosh(fEtaHad), fMassTOFHad)') \
      .Define('fClusterSizeCosLamHe3', 'averageClusterSize(fItsClusterSizeHe3) / cosh(fEtaHe3)') \
      .Define('fClusterSizeCosLamHad', 'averageClusterSize(fItsClusterSizeHad) / cosh(fEtaHad)') \
      .Define('fExpectedClusterSizeHe3', 'ExpectedClusterSizeCosLambdaHe(fPtHe3 * std::cosh(fEtaHe3))') \
      .Define('fExpectedClusterSizeHad', 'ExpectedClusterSizeCosLambdaPr(fPtHad * std::cosh(fEtaHad))') \
      .Define('fNSigmaITSHe3', 'NsigmaITSHe(fPtHe3 * std::cosh(fEtaHe3), fClusterSizeCosLamHe3)') \
      .Define('fNSigmaITSHad', 'NsigmaITSPr(fPtHad * std::cosh(fEtaHad), fClusterSizeCosLamHad)') \
      .Filter(base_selection).Filter(selection) \
      .Define('fKstar', f'Kstar(fPtHe3, fEtaHe3, fPhiHe3, {ParticleMasses["He"]}, fPtHad, fEtaHad, fPhiHad, {ParticleMasses["Pr"]})') \
      .Define('fPxLi', 'fPtHe3 * std::cos(fPhiHe3) + fPtHad * std::cos(fPhiHad)') \
      .Define('fPyLi', 'fPtHe3 * std::sin(fPhiHe3) + fPtHad * std::sin(fPhiHad)') \
      .Define('fPzLi', 'fPtHe3 * std::sinh(fEtaHe3) + fPtHad * std::sinh(fEtaHad)') \
      .Define('fELi', 'fEHe3 + fEHad') \
      .Define('fPLi', 'sqrt(fPxLi*fPxLi + fPyLi*fPyLi + fPzLi*fPzLi)') \
      .Define('fPtLi', 'sqrt(fPxLi*fPxLi + fPyLi*fPyLi)') \
      .Define('fKt', 'fPtLi / 2') \
      .Define('fMt', f'std::sqrt((fPtLi/4)*(fPtLi/4) + {ParticleMasses["Pr"]}*{ParticleMasses["Pr"]})') \
      .Define('fEtaLi', 'std::acosh(fPLi / fELi)') \
      .Define('fPhiLi', 'std::atan2(fPyLi, fPxLi)') \
      .Define('fSignedPtLi', 'fPtLi * fSignHe3') \
      .Define('fMassInvLi', 'std::sqrt(fELi*fELi - fPLi*fPLi)') \
      .Define('fMassTLi', 'std::sqrt(fELi*fELi - fPtLi*fPtLi)') \
      #.Redefine('fNSigmaTPCHe3', 'NsigmaTpcHe(std::abs(fInnerParamTPCHe3), fSignalTPCHe3)') \ after innerparam_tpc is defined

print(f'Selections done!')

histogram_registry = HistogramRegistry()

register_qa_histograms(histogram_registry)
print(f'QA Histograms created!')

register_kstar_histograms(histogram_registry)
print(f'(anti)matter histograms created!')

# Matter
register_kstar_matter_histograms(histogram_registry)
print(f'matter histograms created!')

# Antimatter
register_kstar_antimatter_histograms(histogram_registry)
print(f'antimatter histograms created!')

output_file = config['output_file']
outFile = ROOT.TFile(output_file, "RECREATE")

histogram_registry.prepare_directories(outFile)
histogram_registry.draw_histogram(rdf)
histogram_registry.save_histograms(outFile)

outFile.Close()