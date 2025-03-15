import numpy as np
from ROOT import TDirectory, TFile

from torchic import Dataset, AxisSpec
from torchic.physics.ITS import average_cluster_size

def compute_efficiency(data: Dataset, output_file: TDirectory):

    data['fPtMC'] = data['fPtMC'] / 2.
    rec_data = data.query('fPt > -990', inplace=False)
    gen_data = data

    # do it in two steps for matter and antimatter
    conditions_dict = {'matter': 'fPtMC > 0', 
                       'antimatter': 'fPtMC < 0'}

    axis_spec_pt_rec = AxisSpec(100, 0, 10, 'pt_rec', '#it{p}_{T}^{gen} (GeV/#it{c})')
    axis_spec_pt_gen = AxisSpec(100, 0, 10, 'pt_gen', '#it{p}_{T}^{gen} (GeV/#it{c})')

    for key, condition in conditions_dict.items():
        tmp_rec_data = rec_data.query(condition, inplace=False)
        tmp_gen_data = gen_data.query(condition, inplace=False)

        axis_spec_pt_rec = AxisSpec(100, -10, 0, 'pt_rec', '#it{p}_{T}^{gen} (GeV/#it{c})')
        axis_spec_pt_gen = AxisSpec(100, -10, 0, 'pt_gen', '#it{p}_{T}^{gen} (GeV/#it{c})')
        if key == 'matter':
            axis_spec_pt_rec = AxisSpec(100, 0, 10, 'pt_rec', '#it{p}_{T}^{gen} (GeV/#it{c})')
            axis_spec_pt_gen = AxisSpec(100, 0, 10, 'pt_gen', '#it{p}_{T}^{gen} (GeV/#it{c})')

        h_rec = tmp_rec_data.build_hist('fPtMC', axis_spec_pt_rec)
        h_gen = tmp_gen_data.build_hist('fPtMC', axis_spec_pt_gen)
        h_efficiency = tmp_rec_data.Clone()
        h_efficiency.Reset()

        for ibin in range(1, h_rec.GetNbinsX() + 1):
            rec = h_rec.GetBinContent(ibin)
            gen = h_gen.GetBinContent(ibin)
            eff, eff_err = 0, 0
            if gen > 0:
                eff = rec / gen
                eff_err = np.sqrt(eff * (1 - eff) / gen)
            h_efficiency.SetBinContent(ibin, eff) 
            h_efficiency.SetBinError(ibin, eff_err)

        output_file.cd()
        h_efficiency.Write(f'hEfficiency_{key}')

def data_visual(data: Dataset, output_file: TDirectory):

    axis_spec_pt_rec = AxisSpec(100, -10, 10, 'h_PtRec', '#it{p}_{T}^{rec} (GeV/#it{c})')
    axis_spec_pt_gen = AxisSpec(100, -10, 10, 'h_PtGen', '#it{p}_{T}^{gen} (GeV/#it{c})')
    axis_spec_its = AxisSpec(62, -0.5, 15.5, 'h2_ItsClusterSize', ';#beta#gamma; #LT ITS cluster size #GT #times #LT cos#lambda #GT')

    rec_data = data.query('fPt > -990', inplace=False)
    rec_data['fITSAvgClSize'], __ = average_cluster_size(rec_data['fItsClusterSize'])
    rec_data['fITSClSizeCosLam'] = rec_data['fITSAvgClSize'] / np.cosh(rec_data['fEta'])
    
    h_pt_gen = data.build_hist('fPtMC', axis_spec_pt_gen)
    h_pt_rec = rec_data.build_hist('fPt', axis_spec_pt_rec)
    h2_its_pt = rec_data.build_hist('fPt', 'fItsClSizeCosLam', axis_spec_pt_rec, axis_spec_its)

    output_file.cd()
    h_pt_gen.Write()
    h_pt_rec.Write()
    h2_its_pt.Write()
    

if __name__ == '__main__':

    input_files = {'He': ['/home/galucia/antilithium4/task/MCWorkflowFindables/output/MC_efficiency_he.root'],
                   'Pr': ['/home/galucia/antilithium4/task/MCWorkflowFindables/output/MC_efficiency_pr.root']
                   }

    output_file_path = '/home/galucia/antilithium4/analysis/output/MC/mc.root'
    output_file = TFile(output_file_path, 'recreate')

    for part in ['He', 'Pr']:

        tree_name = 'o2lithium4findmc'
        dir_prefix = 'DF*'
        data = Dataset(input_files[part], tree_name=tree_name, dir_prefix=dir_prefix)

        tdir = output_file.mkdir(part)

        compute_efficiency(data, tdir)
        data_visual(data, tdir)
    
    output_file.Close()