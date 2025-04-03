'''
    Calibration of the ITS cluster size
'''

from torchic import Dataset, AxisSpec
from torchic.physics.ITS import average_cluster_size, expected_cluster_size, sigma_its
from ROOT import TFile

PART_NAME = {
    'He': 'He3',
    'Pr': 'Had',
}
PART_MASS = {
    'He': 2.80923, # GeV/c^2
    'Pr': 0.938272, # GeV/c^2
}

def fit_its_cluster_size(dataset:Dataset, particle:str, outfile:TFile):
    '''
        Fit the cluster size in the ITS for the given particle
        Args:
            dataset: Dataset
                The dataset containing the data
            particle: str
                The particle to fit the cluster size for
            axis: AxisSpec
                The axis to use for the fit
        Returns:
            fit: Fit
                The fit of the cluster size
    '''

    axis_spec_bg = AxisSpec(50, 0, 5, 'bg', ';;')
    axis_spec_clsize = AxisSpec(30, 0, 15, 'cluster_size_cal', ';#beta#gamma;#LT ITS Cluster Size #GT #times cos #LT #lambda #GT')
    h2 = dataset.build_th2(f'fBetaGamma{PART_NAME[particle]}', f'fClSizeCosLam{PART_NAME[particle]}', axis_spec_bg, axis_spec_clsize)

    outfile.cd()
    h2.Write(f'its_cluster_size')


if __name__ == '__main__':

    infile_path = '/data/galucia/lithium_local/same/LHC24ar_pass1_same.root'
    #infile_path = '/data/galucia/lithium_local/same/LHC23_PbPb_pass4_long_same_lsus.root'
    folder_name = 'DF*'
    tree_name = 'O2he3hadtable'

    dataset = Dataset.from_root(infile_path, folder_name=folder_name, tree_name=tree_name, columns=['fItsClusterSizeHe3', 'fItsClusterSizeHad', 'fPtHe3', 'fPtHad', 'fEtaHe3', 'fEtaHad'])
    dataset.eval(f'fBetaGammaHe3 = abs(fPtHe3) * cosh(fEtaHe3) / {PART_MASS["He"]}', inplace=True)
    dataset.eval(f'fBetaGammaHad = abs(fPtHad) * cosh(fEtaHad) / {PART_MASS["Pr"]}', inplace=True)
    dataset['fClSizeITSMeanHe3'], __ = average_cluster_size(dataset['fItsClusterSizeHe3'])
    dataset['fClSizeITSMeanHad'], __ = average_cluster_size(dataset['fItsClusterSizeHad'])
    dataset.eval('fClSizeCosLamHe3 = fClSizeITSMeanHe3 / cosh(fEtaHe3)', inplace=True)
    dataset.eval('fClSizeCosLamHad = fClSizeITSMeanHad / cosh(fEtaHad)', inplace=True)
    
    # sanity check - apply the cut on He3
    dataset['fExpClSizeCosLamHe3'] = expected_cluster_size(dataset['fBetaGammaHe3'], particle='He')
    dataset['fSigmaITSHe3'] = sigma_its(dataset['fClSizeCosLamHe3'], particle='He')
    dataset.eval('fNSigmaITSHe3 = (fClSizeCosLamHe3 - fExpClSizeCosLamHe3) / fSigmaITSHe3', inplace=True)
    dataset.query('fNSigmaITSHe3 > -1.5', inplace=True)

    outfile = TFile.Open('output/ITS_calibration_24.root', 'RECREATE')
    for particle in ['He', 'Pr']:
        outdir = outfile.mkdir(f'{particle}')
        fit_its_cluster_size(dataset, particle, outdir)

    outfile.Close()

