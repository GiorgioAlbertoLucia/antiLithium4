import argparse
import yaml
import numpy as np
import ROOT
from ROOT import TFile, TChain, TList, TTree

from utils.particles import ParticleMasses

ROOT.gROOT.LoadMacro('Common.h++')
from ROOT import NsigmaTpcHe, NsigmaITSHe, NsigmaITSPr, NsigmaTOFPr, averageClusterSize, CorrectPidTrkHe, \
  Kstar, ExpectedClusterSizeCosLambdaHe, ExpectedClusterSizeCosLambdaPr

ROOT.EnableImplicitMT(30)
ROOT.gROOT.SetBatch(True)

base_selection = '(fSignedPtHe3 > 0 && fSignedPtPr > 0) || (fSignedPtHe3 < 0 && fSignedPtPr < 0)'
selections = [
  '((fDeltaPhi/0.006)*(fDeltaPhi/0.006) + (fDeltaEta/0.012)*(fDeltaEta/0.012) > 1)',

  '(std::abs(fEtaHe3) < 0.9)',
  '((fPtHe3 < 2.5) || (fPIDtrkHe3 == 7))',
  #'((0.5 < fChi2TPCHe3) && (fChi2TPCHe3 < 4))',
  '((-1.5 < fNSigmaTPCHe3) && (fNSigmaTPCHe3 < 2.5))',
  #'((-0.9 < fNSigmaTPCHe3) && (fNSigmaTPCHe3 < 1.9))',
  '(std::abs(fDCAxyHe3) < 0.1)',
  '(std::abs(fDCAzHe3) < 1.0)',
  '(fClusterSizeCosLamHe3 > 4)',
  #'(fNSigmaITSHe3 > -1.5)',

  '(std::abs(fEtaPr) < 0.9)',
  #'((0.5 < fChi2TPCPr) && (fChi2TPCPr < 4))',
  '(std::abs(fNSigmaTPCPr) < 2)',
  '((std::abs(fPtPr) < 0.8) || (std::abs(fNSigmaTOFPr) < 2))',
  '(fNSigmaITSPr > -1.5)',
]
#selections = conf['selections']
selection = selections[0]
for sel in selections[1:]:
    selection += (' && ' + sel)
print(f'Selection: {selection}')

############################################################################################################

input_data = ['/data/galucia/lithium_local/same/LHC23_same_chi2tpc.root']
file_data_list = input_data if isinstance(input_data, list) else [input_data]

chainData = TChain("O2lithium4table")
tree_name = "O2lithium4table"

mode = 'DF'  # DataFrame mode

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
      .Define('fSignedPtPr', 'fPtPr') \
      .Define('fSignHe3', 'fPtHe3/std::abs(fPtHe3)') \
      .Redefine('fPtHe3', 'std::abs(fPtHe3)') \
      .Redefine('fPtPr', 'std::abs(fPtPr)') \
      .Redefine('fPtHe3', '(fPIDtrkHe3 == 7) || (fPtHe3 > 2.5) ? fPtHe3 : CorrectPidTrkHe(fPtHe3)') \
      .Define('fSignedPtHe3', 'fPtHe3 * fSignHe3') \
      .Define(f'fEHe3', f'std::sqrt((fPtHe3 * std::cosh(fEtaHe3))*(fPtHe3 * std::cosh(fEtaHe3)) + {ParticleMasses["He"]}*{ParticleMasses["He"]})') \
      .Define(f'fEPr', f'std::sqrt((fPtPr * std::cosh(fEtaPr))*(fPtPr * std::cosh(fEtaPr)) + {ParticleMasses["Pr"]}*{ParticleMasses["Pr"]})') \
      .Define('fDeltaEta', 'fEtaHe3 - fEtaPr') \
      .Define('fDeltaPhi', 'fPhiHe3 - fPhiPr') \
      .Redefine('fInnerParamTPCHe3', 'fInnerParamTPCHe3 * 2') \
      .Redefine('fNSigmaTPCHe3', 'NsigmaTpcHe(std::abs(fInnerParamTPCHe3), fSignalTPCHe3)') \
      .Define('fNSigmaTOFPr', 'NsigmaTOFPr(fPtPr * std::cosh(fEtaPr), fMassTOFPr)') \
      .Define('fClusterSizeCosLamHe3', 'averageClusterSize(fItsClusterSizeHe3) / cosh(fEtaHe3)') \
      .Define('fClusterSizeCosLamPr', 'averageClusterSize(fItsClusterSizePr) / cosh(fEtaPr)') \
      .Define('fExpectedClusterSizeHe3', 'ExpectedClusterSizeCosLambdaHe(fPtHe3)') \
      .Define('fExpectedClusterSizePr', 'ExpectedClusterSizeCosLambdaPr(fPtPr)') \
      .Define('fNSigmaITSHe3', 'NsigmaITSHe(fPtHe3 * std::cosh(fEtaHe3), fClusterSizeCosLamHe3)') \
      .Define('fNSigmaITSPr', 'NsigmaITSPr(fPtPr * std::cosh(fEtaPr), fClusterSizeCosLamPr)') \
      .Filter(base_selection).Filter(selection) 

kstar_min, kstar_max, kstar_bins = 0, 0.8, 40

print(f'Selections done!')
hPhiChi2He3 = rdf.Histo2D(("hPhiChi2He3", ";#phi (rad);#chi^{2}_{TPC} / N_{clusters TPC}", 100, -3.14, 3.14, 100, 0, 10), "fPhiHe3", "fChi2TPCHe3")
hPtChi2He3 = rdf.Histo2D(("hPhiChi2He3", ";#it{p}_{T} (GeV/#it{c});#chi^{2}_{TPC} / N_{clusters TPC}", 100, 0, 10, 100, 0, 10), "fPtHe3", "fChi2TPCHe3")
hPhiChi2Pr = rdf.Histo2D(("hPhiChi2Pr", ";#phi (rad);#chi^{2}_{TPC} / N_{clusters TPC}", 100, -3.14, 3.14, 100, 0, 10), "fPhiPr", "fChi2TPCPr")
print(f'QA Histograms created!')


output_file = 'output/chi2.root'
outFile = ROOT.TFile(output_file, "RECREATE")

qa_dir = outFile.mkdir("QA")
qa_dir.cd()
hPhiChi2He3.Write()
hPtChi2He3.Write()
hPhiChi2Pr.Write()

outFile.Close()