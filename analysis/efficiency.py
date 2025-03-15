from torchic import Dataset, AxisSpec
from torchic.core.histogram import build_efficiency
from src.preprocessing import DataPreprocessor
from copy import deepcopy
from ROOT import TFile

if __name__ == '__main__':

    infile = '/data/galucia/lithium_local/MC/LHC25a4.root'
    dataset = Dataset.from_root(infile, tree_name='O2he3hadtable', folder_name='DF*')
    dataset = dataset.concat(Dataset.from_root(infile, tree_name='O2he3hadtablemc', folder_name='DF*'), axis=1)
    print(f'{dataset.columns=}')

    axis_spec_pt_rec = AxisSpec(100, 0, 10, 'h_PtRec', '#it{p}_{T}^{rec} (GeV/#it{c})')
    axis_spec_pt_gen = AxisSpec(100, 0, 10, 'h_PtGen', '#it{p}_{T}^{gen} (GeV/#it{c})')

    dataset_gen = deepcopy(dataset)
    preprocessor = DataPreprocessor(dataset)
    preprocessor.general_selections()
    preprocessor.define_variables()
    preprocessor.selections_He3()
    preprocessor.selections_Pr()
    dataset_rec = preprocessor.dataset
    dataset_rec.query(f'fPtHe3 > -900', inplace=True)
    dataset_rec.eval('fSignedPtLi = fPtLi * fSignHe3', inplace=True)

    output_file = TFile.Open('/home/galucia/antiLithium4/analysis/output/MC/efficiency.root', 'recreate')

    for type in ['matter', 'antimatter']:
        if type == 'antimatter':
            dataset_gen[f'fSignedPtMC'] = -dataset_gen[f'fSignedPtMC']
            dataset_rec[f'fSignedPtLi'] = -dataset_rec[f'fSignedPtLi']

        h_pt_gen = dataset_gen.build_th1(f'fSignedPtMC', axis_spec_pt_gen)
        h_pt_rec = dataset_rec.build_th1(f'fSignedPtLi', axis_spec_pt_rec)
        h_efficiency = build_efficiency(h_pt_gen, h_pt_rec, f'efficiency_{type}', '#it{p}_{T}^{rec} (GeV/#it{c})', 'Efficiency')

        output_file.cd()
        h_efficiency.Write()
        del h_pt_gen, h_pt_rec, h_efficiency
    
    output_file.Close()