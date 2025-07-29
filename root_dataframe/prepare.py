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

#parser = argparse.ArgumentParser(description='Configure the parameters of the script.')
#parser.add_argument('--config-file', dest='config_file', help='path to the YAML file with configuration.', default='')
#args = parser.parse_args()
#if args.config_file == '':
#    print('** No config file provided. Exiting. **')
#    exit()

#with open(args.config_file, 'r') as stream:
#  confFile = yaml.safe_load(stream)
#  conf = confFile["data_preparation"]

base_selection = '(fSignedPtHe3 > 0 && fSignedPtHad > 0) || (fSignedPtHe3 < 0 && fSignedPtHad < 0)'
selections = [
  '((fDeltaPhi/0.006)*(fDeltaPhi/0.006) + (fDeltaEta/0.012)*(fDeltaEta/0.012) > 1)',

  '(std::abs(fEtaHe3) < 0.9)',
  '((fPtHe3 < 2.5) || (fPIDtrkHe3 == 7))',
  '((0.5 < fChi2TPCHe3) && (fChi2TPCHe3 < 4))',
  '((-1.5 < fNSigmaTPCHe3) && (fNSigmaTPCHe3 < 2.5))',
  #'((-0.9 < fNSigmaTPCHe3) && (fNSigmaTPCHe3 < 1.9))',
  '(std::abs(fDCAxyHe3) < 0.1)',
  '(std::abs(fDCAzHe3) < 1.0)',
  '(fClusterSizeCosLamHe3 > 4)',
  '(fNSigmaITSHe3 > -1.5)',

  '(std::abs(fEtaHad) < 0.9)',
  '((0.5 < fChi2TPCHad) && (fChi2TPCHad < 4))',
  '(std::abs(fNSigmaTPCHad) < 2)',
  '((std::abs(fPtHad) < 0.8) || (std::abs(fNSigmaTOFHad) < 2))',
  '(fNSigmaITSHad > -1.5)',
]
#selections = conf['selections']
selection = selections[0]
for sel in selections[1:]:
    selection += (' && ' + sel)
print(f'Selection: {selection}')

############################################################################################################

#input_data = ['input/mc.root',
#              ]
input_data = ['input/LHC23_PbPb_pass4_long_same_lsus_merged.root',
              'input/LHC24ar_pass1_same_merged.root',
              'input/LHC24as_pass1_same_merged.root'
              ]
#input_data = ['/data/galucia/lithium_local/mixing/LHC23_PbPb_pass4_long_mixing_lsus.root',
#              '/data/galucia/lithium_local/mixing/LHC24ar_pass1_mixing_lsus.root',
#              '/data/galucia/lithium_local/mixing/LHC24as_pass1_mixing_lsus.root',
#              ]
file_data_list = input_data if isinstance(input_data, list) else [input_data]

chainData = TChain("O2he3hadtable")
tree_name = "O2he3hadtable"
#tree_name = "MixedTree"

mode = 'DF'  # DataFrame mode
#mode = 'tree'

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
      .Redefine('fNSigmaTPCHe3', 'NsigmaTpcHe(std::abs(fInnerParamTPCHe3), fSignalTPCHe3)') \
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

kstar_min, kstar_max, kstar_bins = 0, 0.8, 80

print(f'Selections done!')
pt_bins = np.linspace(-10, 10, 200)
h2PtNSigmaTPCHe = rdf.Histo2D(("h2PtNSigmaTPCHe", ";#it{p}_{T} (GeV/#it{c});n#sigma_{TPC}", len(pt_bins) - 1, pt_bins, 100, -4, 4), "fSignedPtHe3", "fNSigmaTPCHe3")
h2PtNSigmaTPCPr = rdf.Histo2D(("h2PtNSigmaTPCPr", ";#it{p}_{T} (GeV/#it{c});n#sigma_{TPC}", len(pt_bins) - 1, pt_bins, 100, -4, 4), "fSignedPtHad", "fNSigmaTPCHad")
h2PtClusterSizeHe = rdf.Histo2D(("h2PtClusterSizeHe", ";#it{p}_{T} (GeV/#it{c});Cluster size", len(pt_bins) - 1, pt_bins, 90, 0, 15), "fSignedPtHe3", "fClusterSizeCosLamHe3")
h2PtExpectedClusterSizeHe = rdf.Histo2D(("h2PtExpectedClusterSizeHe", ";#it{p}_{T} (GeV/#it{c});Expected cluster size", len(pt_bins) - 1, pt_bins, 90, 0, 15), "fSignedPtHe3", "fExpectedClusterSizeHe3")
h2PtClusterSizePr = rdf.Histo2D(("h2PtClusterSizePr", ";#it{p}_{T} (GeV/#it{c});Cluster size", len(pt_bins) - 1, pt_bins, 90, 0, 15), "fSignedPtHad", "fClusterSizeCosLamHad")
h2PtExpectedClusterSizePr = rdf.Histo2D(("h2PtExpectedClusterSizePr", ";#it{p}_{T} (GeV/#it{c});Expected cluster size", len(pt_bins) - 1, pt_bins, 90, 0, 15), "fSignedPtHad", "fExpectedClusterSizeHad")
h2PtNSigmaITSHe = rdf.Histo2D(("h2PtNSigmaITSHe", ";#it{p}_{T} (GeV/#it{c});n#sigma_{ITS}", len(pt_bins) - 1, pt_bins, 100, -4, 4), "fSignedPtHe3", "fNSigmaITSHe3")
h2PtNSigmaITSPr = rdf.Histo2D(("h2PtNSigmaITSPr", ";#it{p}_{T} (GeV/#it{c});n#sigma_{ITS}", len(pt_bins) - 1, pt_bins, 100, -4, 4), "fSignedPtHad", "fNSigmaITSHad")
h2DeltaEtaDeltaPhi = rdf.Histo2D(("h2DeltaEtaDeltaPhi", ";#Delta#eta;#Delta#phi (rad)", 100, -0.1, 0.1, 100, -0.1, 0.1), "fDeltaEta", "fDeltaPhi")
hEta = rdf.Histo1D(("hEtaHe", ";#it{k}^{*} (GeV/#it{c});", 100, -1, 1.), "fEtaHe3")
hKstar = rdf.Histo1D(("hKstar", ";#it{k}^{*} (GeV/#it{c});", kstar_bins, kstar_min, kstar_max), "fKstar")
hCentrality = rdf.Histo1D(("hCentrality", ";Centrality (ft0c);", 100, 0, 100), "fCentralityFT0C")
hInvariantMass = rdf.Histo1D(("hInvariantMass", ";Invariant mass (GeV/#it{c}^{2});", 400, 3.747, 3.947), "fMassInvLi")
hPtDCAxyHe3 = rdf.Histo2D(("hPtDCAxyHe3", ";#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)", 100, 0, 10, 100, -0.1, 0.1), "fSignedPtHe3", "fDCAxyHe3")
hPtDCAzHe3 = rdf.Histo2D(("hPtDCAzHe3", ";#it{p}_{T} (GeV/#it{c});DCA_{z} (cm)", 100, 0, 10, 100, -1.0, 1.0), "fSignedPtHe3", "fDCAzHe3")
hPtDCAxyHad = rdf.Histo2D(("hPtDCAxyHad", ";#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)", 100, 0, 10, 100, -0.1, 0.1), "fSignedPtHad", "fDCAxyHad")
hPtDCAzHad = rdf.Histo2D(("hPtDCAzHad", ";#it{p}_{T} (GeV/#it{c});DCA_{z} (cm)", 100, 0, 10, 100, -1.0, 1.0), "fSignedPtHad", "fDCAzHad")

hPhiChi2He3 = rdf.Histo2D(("hPhiChi2He3", ";#phi (rad);#chi^{2}_{TPC} / N_{clusters TPC}", 100, -3.14, 3.14, 100, 0, 10), "fPhiHe3", "fChi2TPCHe3")
hPhiChi2Had = rdf.Histo2D(("hPhiChi2Had", ";#phi (rad);#chi^{2}_{TPC} / N_{clusters TPC}", 100, -3.14, 3.14, 100, 0, 10), "fPhiHad", "fChi2TPCHad")
print(f'QA Histograms created!')

hCentralityKstar = rdf.Histo2D(("hCentralityKstar", ";Centrality (ft0c);#it{k}^{*} (GeV/#it{c})", 100, 0, 100, 100, 0, 1.), "fCentralityFT0C", "fKstar")
hKstar010 = rdf.Filter("fCentralityFT0C < 10").Histo1D(("hKstar010", ";#it{k}^{*} (GeV/#it{c});", kstar_bins, kstar_min, kstar_max), "fKstar")
hKstar010Limited = rdf.Filter("1.1 < fKt && fKt < 2 && fCentralityFT0C < 10").Histo1D(("hKstar010Limited", ";#it{k}^{*} (GeV/#it{c});", kstar_bins, kstar_min, kstar_max), "fKstar")
hKt010 = rdf.Filter("fCentralityFT0C < 10").Histo1D(("hKt010", ";#it{k}_{T} (GeV/#it{c});", 500, 0, 5), "fKt")
hKt010Limited = rdf.Filter("1.1 < fKt && fKt < 2 && fCentralityFT0C < 10").Histo1D(("hKt010Limited", ";#it{k}_{T} (GeV/#it{c});", 500, 0, 5), "fKt")
hMt010 = rdf.Filter("fCentralityFT0C < 10").Histo1D(("hMt010", ";#it{m}_{T} (GeV/#it{c});", 500, 0, 5), "fMt")
hKstar1030 = rdf.Filter("10 <= fCentralityFT0C && fCentralityFT0C < 30").Histo1D(("hKstar1030", ";#it{k}^{*} (GeV/#it{c});", kstar_bins, kstar_min, kstar_max), "fKstar")
hKstar1030Limited = rdf.Filter("1.1 < fKt && fKt < 2 && 10 <= fCentralityFT0C && fCentralityFT0C < 30").Histo1D(("hKstar1030Limited", ";#it{k}^{*} (GeV/#it{c});", kstar_bins, kstar_min, kstar_max), "fKstar")
hKt1030 = rdf.Filter("10 <= fCentralityFT0C && fCentralityFT0C < 30").Histo1D(("hKt1030", ";#it{k}_{T} (GeV/#it{c});", 500, 0, 5), "fKt")
hKt1030Limited = rdf.Filter("1.1 < fKt && fKt < 2 && 10 <= fCentralityFT0C && fCentralityFT0C < 30").Histo1D(("hKt1030Limited", ";#it{k}_{T} (GeV/#it{c});", 500, 0, 5), "fKt")
hMt1030 = rdf.Filter("10 <= fCentralityFT0C && fCentralityFT0C < 30").Histo1D(("hMt1030", ";#it{m}_{T} (GeV/#it{c});", 500, 0, 5), "fMt")
hKstar3050 = rdf.Filter("30 <= fCentralityFT0C && fCentralityFT0C < 50").Histo1D(("hKstar3050", ";#it{k}^{*} (GeV/#it{c});", kstar_bins, kstar_min, kstar_max), "fKstar")
hKstar3050Limited = rdf.Filter("1.1 < fKt && fKt < 2 && 30 <= fCentralityFT0C && fCentralityFT0C < 50").Histo1D(("hKstar3050Limited", ";#it{k}^{*} (GeV/#it{c});", kstar_bins, kstar_min, kstar_max), "fKstar")
hKt3050 = rdf.Filter("30 <= fCentralityFT0C && fCentralityFT0C < 50").Histo1D(("hKt3050", ";#it{k}_{T} (GeV/#it{c});", 500, 0, 5), "fKt")
hKt3050Limited = rdf.Filter("1.1 < fKt && fKt < 2 && 30 <= fCentralityFT0C && fCentralityFT0C < 50").Histo1D(("hKt3050Limited", ";#it{k}_{T} (GeV/#it{c});", 500, 0, 5), "fKt")
hMt3050 = rdf.Filter("30 <= fCentralityFT0C && fCentralityFT0C < 50").Histo1D(("hMt3050", ";#it{m}_{T} (GeV/#it{c});", 500, 0, 5), "fMt")
print(f'(anti)matter histograms created!')

# Matter
hCentralityKstarMatter = rdf.Filter("fSignedPtHe3 > 0").Histo2D(("hCentralityKstarMatter", ";Centrality (ft0c);#it{k}^{*} (GeV/#it{c})", 100, 0, 100, 100, 0, 1.), "fCentralityFT0C", "fKstar")
hKstarMatter = rdf.Filter("fSignedPtHe3 > 0").Histo1D(("hKstarMatter", ";#it{k}^{*} (GeV/#it{c});", kstar_bins, kstar_min, kstar_max), "fKstar")
hKstar010Matter = rdf.Filter("fCentralityFT0C < 10 && fSignedPtHe3 > 0").Histo1D(("hKstar010Matter", ";#it{k}^{*} (GeV/#it{c});", kstar_bins, kstar_min, kstar_max), "fKstar")
hKstar010MatterLimited = rdf.Filter("1.1 < fKt && fKt < 2 && fCentralityFT0C < 10 && fSignedPtHe3 > 0").Histo1D(("hKstar010MatterLimited", ";#it{k}^{*} (GeV/#it{c});", kstar_bins, kstar_min, kstar_max), "fKstar")
hKt010Matter = rdf.Filter("fCentralityFT0C < 10 && fSignedPtHe3 > 0").Histo1D(("hKt010Matter", ";#it{k}_{T} (GeV/#it{c});", 500, 0, 5), "fKt")
hKt010MatterLimited = rdf.Filter("1.1 < fKt && fKt < 2 && fCentralityFT0C < 10 && fSignedPtHe3 > 0").Histo1D(("hKt010MatterLimited", ";#it{k}_{T} (GeV/#it{c});", 500, 0, 5), "fKt")
hMt010Matter = rdf.Filter("fCentralityFT0C < 10 && fSignedPtHe3 > 0").Histo1D(("hMt010Matter", ";#it{m}_{T} (GeV/#it{c});", 500, 0, 5), "fMt")
hKstar1030Matter = rdf.Filter("10 <= fCentralityFT0C && fCentralityFT0C < 30 && fSignedPtHe3 > 0").Histo1D(("hKstar1030Matter", ";#it{k}^{*} (GeV/#it{c});", kstar_bins, kstar_min, kstar_max), "fKstar")
hKstar1030MatterLimited = rdf.Filter("1.1 < fKt && fKt < 2 && 10 <= fCentralityFT0C && fCentralityFT0C < 30 && fSignedPtHe3 > 0").Histo1D(("hKstar1030MatterLimited", ";#it{k}^{*} (GeV/#it{c});", kstar_bins, kstar_min, kstar_max), "fKstar")
hKt1030Matter = rdf.Filter("10 <= fCentralityFT0C && fCentralityFT0C < 30 && fSignedPtHe3 > 0").Histo1D(("hKt1030Matter", ";#it{k}_{T} (GeV/#it{c});", 500, 0, 5), "fKt")
hKt1030MatterLimited = rdf.Filter("1.1 < fKt && fKt < 2 && 10 <= fCentralityFT0C && fCentralityFT0C < 30 && fSignedPtHe3 > 0").Histo1D(("hKt1030MatterLimited", ";#it{k}_{T} (GeV/#it{c});", 500, 0, 5), "fKt")
hMt1030Matter = rdf.Filter("10 <= fCentralityFT0C && fCentralityFT0C < 30 && fSignedPtHe3 > 0").Histo1D(("hMt1030Matter", ";#it{m}_{T} (GeV/#it{c});", 500, 0, 5), "fMt")
hKstar3050Matter = rdf.Filter("30 <= fCentralityFT0C && fCentralityFT0C < 50 && fSignedPtHe3 > 0").Histo1D(("hKstar3050Matter", ";#it{k}^{*} (GeV/#it{c});", kstar_bins, kstar_min, kstar_max), "fKstar")
hKstar3050MatterLimited = rdf.Filter("1.1 < fKt && fKt < 2 && 30 <= fCentralityFT0C && fCentralityFT0C < 50 && fSignedPtHe3 > 0").Histo1D(("hKstar3050MatterLimited", ";#it{k}^{*} (GeV/#it{c});", kstar_bins, kstar_min, kstar_max), "fKstar")
hKt3050Matter = rdf.Filter("30 <= fCentralityFT0C && fCentralityFT0C < 50 && fSignedPtHe3 > 0").Histo1D(("hKt3050Matter", ";#it{k}_{T} (GeV/#it{c});", 500, 0, 5), "fKt")
hKt3050MatterLimited = rdf.Filter("1.1 < fKt && fKt < 2 && 30 <= fCentralityFT0C && fCentralityFT0C < 50 && fSignedPtHe3 > 0").Histo1D(("hKt3050MatterLimited", ";#it{k}_{T} (GeV/#it{c});", 500, 0, 5), "fKt")
hMt3050Matter = rdf.Filter("30 <= fCentralityFT0C && fCentralityFT0C < 50 && fSignedPtHe3 > 0").Histo1D(("hMt3050Matter", ";#it{m}_{T} (GeV/#it{c});", 500, 0, 5), "fMt")
hInvariantMassMatter = rdf.Filter("fSignedPtHe3 > 0").Histo1D(("hInvariantMassMatter", ";Invariant mass (GeV/#it{c}^{2});", 400, 3.747, 3.947), "fMassInvLi")
print(f'matter histograms created!')

# Antimatter
hCentralityKstarAntimatter = rdf.Filter("fSignedPtHe3 < 0").Histo2D(("hCentralityKstarAntimatter", ";Centrality (ft0c);#it{k}^{*} (GeV/#it{c})", 100, 0, 100, 100, 0, 1.), "fCentralityFT0C", "fKstar")
hKstarAntimatter = rdf.Filter("fSignedPtHe3 < 0").Histo1D(("hKstarAntimatter", ";#it{k}^{*} (GeV/#it{c});", kstar_bins, kstar_min, kstar_max), "fKstar")
hKstar010Antimatter = rdf.Filter("fCentralityFT0C < 10 && fSignedPtHe3 < 0").Histo1D(("hKstar010Antimatter", ";#it{k}^{*} (GeV/#it{c});", kstar_bins, kstar_min, kstar_max), "fKstar")
hKstar010AntimatterLimited = rdf.Filter("1.1 < fKt && fKt < 2 && fCentralityFT0C < 10 && fSignedPtHe3 < 0").Histo1D(("hKstar010AntimatterLimited", ";#it{k}^{*} (GeV/#it{c});", kstar_bins, kstar_min, kstar_max), "fKstar")
hKt010Antimatter = rdf.Filter("fCentralityFT0C < 10 && fSignedPtHe3 < 0").Histo1D(("hKt010Antimatter", ";#it{k}_{T} (GeV/#it{c});", 500, 0, 5), "fKt")
hKt010AntimatterLimited = rdf.Filter("1.1 < fKt && fKt < 2 && fCentralityFT0C < 10 && fSignedPtHe3 < 0").Histo1D(("hKt010AntimatterLimited", ";#it{k}_{T} (GeV/#it{c});", 500, 0, 5), "fKt")
hMt010Antimatter = rdf.Filter("fCentralityFT0C < 10 && fSignedPtHe3 < 0").Histo1D(("hMt010Antimatter", ";#it{m}_{T} (GeV/#it{c});", 500, 0, 5), "fMt")
hKstar1030Antimatter = rdf.Filter("10 <= fCentralityFT0C && fCentralityFT0C < 30 && fSignedPtHe3 < 0").Histo1D(("hKstar1030Antimatter", ";#it{k}^{*} (GeV/#it{c});", kstar_bins, kstar_min, kstar_max), "fKstar")
hKstar1030AntimatterLimited = rdf.Filter("1.1 < fKt && fKt < 2 && 10 <= fCentralityFT0C && fCentralityFT0C < 30 && fSignedPtHe3 < 0").Histo1D(("hKstar1030AntimatterLimited", ";#it{k}^{*} (GeV/#it{c});", kstar_bins, kstar_min, kstar_max), "fKstar")
hKt1030Antimatter = rdf.Filter("10 <= fCentralityFT0C && fCentralityFT0C < 30 && fSignedPtHe3 < 0").Histo1D(("hKt1030Antimatter", ";#it{k}_{T} (GeV/#it{c});", 500, 0, 5), "fKt")
hKt1030AntimatterLimited = rdf.Filter("1.1 < fKt && fKt < 2 && 10 <= fCentralityFT0C && fCentralityFT0C < 30 && fSignedPtHe3 < 0").Histo1D(("hKt1030AntimatterLimited", ";#it{k}_{T} (GeV/#it{c});", 500, 0, 5), "fKt")
hMt1030Antimatter = rdf.Filter("10 <= fCentralityFT0C && fCentralityFT0C < 30 && fSignedPtHe3 < 0").Histo1D(("hMt1030Antimatter", ";#it{m}_{T} (GeV/#it{c});", 500, 0, 5), "fMt")
hKstar3050Antimatter = rdf.Filter("30 <= fCentralityFT0C && fCentralityFT0C < 50 && fSignedPtHe3 < 0").Histo1D(("hKstar3050Antimatter", ";#it{k}^{*} (GeV/#it{c});", kstar_bins, kstar_min, kstar_max), "fKstar")
hKstar3050AntimatterLimited = rdf.Filter("1.1 < fKt && fKt < 2 && 30 <= fCentralityFT0C && fCentralityFT0C < 50 && fSignedPtHe3 < 0").Histo1D(("hKstar3050AntimatterLimited", ";#it{k}^{*} (GeV/#it{c});", kstar_bins, kstar_min, kstar_max), "fKstar")
hKt3050Antimatter = rdf.Filter("30 <= fCentralityFT0C && fCentralityFT0C < 50 && fSignedPtHe3 < 0").Histo1D(("hKt3050Antimatter", ";#it{k}_{T} (GeV/#it{c});", 500, 0, 5), "fKt")
hKt3050AntimatterLimited = rdf.Filter("1.1 < fKt && fKt < 2 && 30 <= fCentralityFT0C && fCentralityFT0C < 50 && fSignedPtHe3 < 0").Histo1D(("hKt3050AntimatterLimited", ";#it{k}_{T} (GeV/#it{c});", 500, 0, 5), "fKt")
hMt3050Antimatter = rdf.Filter("30 <= fCentralityFT0C && fCentralityFT0C < 50 && fSignedPtHe3 < 0").Histo1D(("hMt3050Antimatter", ";#it{m}_{T} (GeV/#it{c});", 500, 0, 5), "fMt")
hInvariantMassAntimatter = rdf.Filter("fSignedPtHe3 < 0").Histo1D(("hInvariantMassAntimatter", ";Invariant mass (GeV/#it{c}^{2});", 400, 3.747, 3.947), "fMassInvLi")
print(f'antimatter histograms created!')

#hDCAxyDCaz = rdf.Histo2D(("hDCAxyDCaz", ";DCA_{xy} (cm);DCA_{z} (cm)", n_bins_dca, -max_abs_dca, max_abs_dca, n_bins_dca, -max_abs_dca, max_abs_dca), "fDCAxy", "fDCAz")

output_file = 'output/same_event.root'
#output_file = 'output/mixed_event.root'
#output_file = 'output/mc.root'
outFile = ROOT.TFile(output_file, "RECREATE")

qa_dir = outFile.mkdir("QA")
qa_dir.cd()
hEta.Write()
h2PtNSigmaTPCHe.Write()
h2PtNSigmaTPCPr.Write()
h2PtNSigmaITSHe.Write()
h2PtNSigmaITSPr.Write()
h2PtClusterSizeHe.Write()
h2PtClusterSizePr.Write()
h2PtExpectedClusterSizeHe.Write()
h2PtExpectedClusterSizePr.Write()
h2DeltaEtaDeltaPhi.Write()
hPhiChi2He3.Write()
hPhiChi2Had.Write()
hPtDCAxyHe3.Write()
hPtDCAzHe3.Write()
hPtDCAxyHad.Write()
hPtDCAzHad.Write()

kstar_dir = outFile.mkdir("kstar")
kstar_dir.cd()
hKstar.Write()
hCentrality.Write()
hCentralityKstar.Write()
hKstar010.Write()
hKstar1030.Write()
hKstar3050.Write()
hKt010.Write()
hKt1030.Write()
hKt3050.Write()
hKt010Limited.Write()
hKt1030Limited.Write()
hKt3050Limited.Write()
hKstar010Limited.Write()
hKt010Limited.Write()
hKstar1030Limited.Write()
hKt1030Limited.Write()
hKstar3050Limited.Write()
hMt010.Write()
hMt1030.Write()
hMt3050.Write()

invmass_dir = outFile.mkdir("InvariantMass")
invmass_dir.cd()
hInvariantMass.Write()


kstar_dir_matter = outFile.mkdir("kstarMatter")
kstar_dir_matter.cd()
hCentralityKstarMatter.Write()
hKstarMatter.Write()
hKstar010Matter.Write()
hKstar1030Matter.Write()
hKstar3050Matter.Write()
hKt010Matter.Write()
hKt1030Matter.Write()
hKt3050Matter.Write()
hKt010MatterLimited.Write()
hKt1030MatterLimited.Write()
hKt3050MatterLimited.Write()
hKstar010MatterLimited.Write()
hKstar1030MatterLimited.Write()
hKstar3050MatterLimited.Write() 
hMt010Matter.Write()
hMt1030Matter.Write()
hMt3050Matter.Write()

inv_mass_dir_matter = outFile.mkdir("InvariantMassMatter")
inv_mass_dir_matter.cd()
hInvariantMassMatter.Write()  # This histogram is not defined in the script, so it is commented out.

kstar_dir_antimatter = outFile.mkdir("kstarAntimatter")
kstar_dir_antimatter.cd()
hCentralityKstarAntimatter.Write()
hKstarAntimatter.Write()
hKstar010Antimatter.Write()
hKstar1030Antimatter.Write()
hKstar3050Antimatter.Write()
hKt010Antimatter.Write()
hKt1030Antimatter.Write()
hKt3050Antimatter.Write()
hKt010AntimatterLimited.Write()
hKt1030AntimatterLimited.Write()
hKt3050AntimatterLimited.Write()
hKstar010AntimatterLimited.Write()
hKstar1030AntimatterLimited.Write()
hKstar3050AntimatterLimited.Write()
hMt010Antimatter.Write()
hMt1030Antimatter.Write()
hMt3050Antimatter.Write()

inv_mass_dir_antimatter = outFile.mkdir("InvariantMassAntimatter")
inv_mass_dir_antimatter.cd()
hInvariantMassAntimatter.Write()  # This histogram is not defined in the script, so it is commented out.

outFile.Close()