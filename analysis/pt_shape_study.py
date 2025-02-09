'''
    Code to run calibration of ITS and TPC parametrisations
'''

from ROOT import TFile
from torchic import Dataset, AxisSpec

if __name__ == '__main__':

    infile_path = '/Users/glucia/Projects/ALICE/data/lithium/same/LHC24as_pass1_same.root'
    folder_name = 'DF*'
    outfile_path = f'/Users/glucia/Projects/ALICE/antiLithium4/analysis/output/LHC24PbPb/pt_shape.root'
    outfile = TFile.Open(outfile_path, 'recreate')

    dataset_table = Dataset.from_root(infile_path, 'O2he3hadtable', folder_name, columns=['fPtHad', 'fPtHe3'])
    dataset_mult = Dataset.from_root(infile_path, 'O2he3hadmult', folder_name, columns=['fCollisionId'])
    dataset = dataset_table.concat(dataset_mult, axis=1)

    collisions = dataset['fCollisionId'].unique()
    n_part = {collision: dataset.query(f'fCollisionId == {collision}', inplace=False).shape[0] for collision in collisions}
    dataset['fNPart'] = dataset['fCollisionId'].apply(lambda x: n_part[x])
    
    axis_spec_pt = AxisSpec(50, -5, 0, 'pt', ';#it{p}_{T} (GeV/c);')
    axis_spec_npart_he3 = AxisSpec(10, 0, 9000, 'npart_He3', ';#it{p}_{T} (GeV/c);N_{part};')
    axis_spec_npart_had = AxisSpec(10, 0, 9000, 'npart_Had', ';#it{p}_{T} (GeV/c);N_{part};')
    h2_pt_he3 = dataset.build_th2('fPtHe3', 'fNPart', axis_spec_pt, axis_spec_npart_he3)
    h2_pt_had = dataset.build_th2('fPtHad', 'fNPart', axis_spec_pt, axis_spec_npart_had)

    outfile.cd()
    h2_pt_he3.Write()
    h2_pt_had.Write()
    outfile.Close()
