'''
    Preliminar analysis of the montecarlo to check the variables for the machine learning algorithms
'''

import sys
sys.path.append('..')
#from analysis.src.preprocessing import Preprocessor

import uproot
from ROOT import TFile, TCanvas, gStyle, TPaveText
from torchic import Dataset, AxisSpec

def load_data() -> Dataset:

    infile = '/home/galucia/antiLithium4/task/MCWorkflowAnalysis/AO2D_lit_mc.root'
    tree_names = ['O2he3hadtable', 'O2he3hadtablemc', 'O2he3hadmult', 'O2he3hadmcext']
    datasets = []

    for tree_name in tree_names:
        _dataset = Dataset.from_root(infile, tree_name=tree_name, folder_name='DF*')
        datasets.append(_dataset)

    dataset = datasets[0]
    for _dataset in datasets[1:]:
        dataset = dataset.concat(_dataset, axis=1)

    return dataset

def plot(dataset: Dataset, outfile_name: str):

    subsets = ['primary_he3', 'primary_p', 'li4', 'secondary_he3', 'secondary_p', 'secondaries']
    dataset.add_subset('primary_he3', dataset['fIsHe3Primary'] == True)
    dataset.add_subset('primary_p', dataset['fIsHadPrimary'] == True)
    dataset.add_subset('li4', dataset['fIsMotherLi4'] == True)
    dataset.add_subset('secondary_he3', dataset['fIsHe3Primary'] == False)
    dataset.add_subset('secondary_p', dataset['fIsHadPrimary'] == False)
    dataset.add_subset('secondaries', (dataset['fIsHadPrimary'] == False) & (dataset['fIsHadPrimary'] == False))
    
    outfile = TFile.Open(outfile_name, 'recreate')
    for subset in subsets:
    
        axis_spec_dcaxy = AxisSpec(400, -0.2, 0.2, 'dcaxy', ';DCAxy;')
        axis_spec_dcaz = AxisSpec(400, -0.2, 0.2, 'dcaz', ';DCAxy;')
        axis_spec_dcapair = AxisSpec(200, 0., 0.2, 'dcapair', ';DCAdau;')

        h_dcaxyhe3 = dataset.build_th1('fDCAxyHe3', axis_spec_dcaxy, subset=subset)
        h_dcaxyhe3.SetName('dcaxyhe3')
        h_dcazhe3 = dataset.build_th1('fDCAzHe3', axis_spec_dcaz, subset=subset)
        h_dcazhe3.SetName('dcazhe3')
        h_dcaxyp = dataset.build_th1('fDCAxyHad', axis_spec_dcaxy, subset=subset)
        h_dcaxyp.SetName('dcaxyp')
        h_dcazp = dataset.build_th1('fDCAzHad', axis_spec_dcaz, subset=subset)
        h_dcazp.SetName('dcazp')
        h_dcapair = dataset.build_th1('fDCApair', axis_spec_dcapair, subset=subset)
        
        outdir = outfile.mkdir(subset)
        outdir.cd()
        h_dcaxyhe3.Write('dcaxyhe3')
        h_dcazhe3.Write('dcazhe3')
        h_dcaxyp.Write('dcaxyp')
        h_dcazp.Write('dcazp')
        h_dcapair.Write('dcapair')

def plot_canvas(outfile_name:str, outfolder_name:str):

    subsets = {
        'he': ['primary_he3', 'li4', 'secondary_he3'],
        'pr': ['primary_p', 'li4', 'secondary_p']
    }
    hist_names = ['dcaxyhe3', 'dcazhe3', 'dcaxyp', 'dcazp', 'dcapair']

    file = TFile.Open(outfile_name, 'read')
    gStyle.SetOptStat(0)

    for hist_name in hist_names:
        subset = None
        if 'he3' in hist_name:
            subset = subsets['he']
        elif 'p' in hist_name:
            subset = subsets['pr']
        else:
            subset = set(subsets['he'] + subsets['pr'])
        
        canvas = TCanvas(f'c_{hist_name}', f'{hist_name}', 800, 600)

        canvas.DrawFrame(-0.2, 1, 0.2, 5e5, f'{hist_name};DCA;Counts')
        text = TPaveText(0.05, 5e2, 0.19, 1e4)
        text.SetFillColor(0)

        for set in subset:
            hist = file.Get(f'{set}/{hist_name}')
            hist.SetLineColor(subset.index(set) + 1)
            hist.SetName(f'{set}_{hist_name}')
            hist.Draw('same')
            text.AddText(f'{set} RMS: {hist.GetRMS():.4f}')
        
        canvas.SetLogy()
        text.Draw('same')
        canvas.BuildLegend(0.7, 0.7, 0.89, 0.89)
        canvas.SaveAs(f'{outfolder_name}/{hist_name}.pdf')



if __name__ == '__main__':

    outfile_name = 'output/analysis.root'
    outfolder_name = 'output'

    #dataset = load_data()
    #print(f'{dataset.columns=}')
    #plot(dataset, outfile_name)
    plot_canvas(outfile_name, outfolder_name)