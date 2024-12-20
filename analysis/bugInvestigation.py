from torchic import Dataset
import sys
sys.path.append('/home/galucia/antiLithium4/analysis')
from utils.particles import ParticleMasses

if __name__ == '__main__':

    dataset = Dataset.from_root('/data/galucia/lithium_local/same/LHC23_PbPb_pass4_long_same_us.root', folder_name='DF*', tree_name='O2he3hadtable')

    # cut in pseudorapidity
    dataset.query('-0.9 < fEtaHe3 < 0.9', inplace=True)
    dataset.query('-0.9 < fEtaHad < 0.9', inplace=True)

    ## definition of reconstructed variables
    dataset.eval('fSignedPtHe3 = fPtHe3', inplace=True)
    dataset.eval('fSignedPtHad = fPtHad', inplace=True)
    dataset.eval('fPtHe3 = abs(fPtHe3)', inplace=True)
    dataset.eval('fPtHad = abs(fPtHad)', inplace=True)
    dataset.eval('fSignHe3 = fSignedPtHe3/fPtHe3', inplace=True)
    dataset.eval('fSignHad = fSignedPtHad/fPtHad', inplace=True)

    dataset.eval('fPxHe3 = fPtHe3 * cos(fPhiHe3)', inplace=True)
    dataset.eval('fPxHad = fPtHad * cos(fPhiHad)', inplace=True)
    dataset.eval('fPyHe3 = fPtHe3 * sin(fPhiHe3)', inplace=True)
    dataset.eval('fPyHad = fPtHad * sin(fPhiHad)', inplace=True)
    dataset.eval('fPzHe3 = fPtHe3 * sinh(fEtaHe3)', inplace=True)
    dataset.eval('fPzHad = fPtHad * sinh(fEtaHad)', inplace=True)
    dataset.eval('fPHe3 = fPtHe3 * cosh(fEtaHe3)', inplace=True)
    dataset.eval('fPHad = fPtHad * cosh(fEtaHad)', inplace=True)
    dataset.eval(f'fEHe3 = sqrt(fPHe3**2 + {ParticleMasses["He"]}**2)', inplace=True)
    dataset.eval(f'fEHad = sqrt(fPHad**2 + {ParticleMasses["Pr"]}**2)', inplace=True)

    # invariant mass 
    dataset.eval('fPxLi = fPxHe3 + fPxHad', inplace=True)
    dataset.eval('fPyLi = fPyHe3 + fPyHad', inplace=True)
    dataset.eval('fPzLi = fPzHe3 + fPzHad', inplace=True)
    dataset.eval('fELi = fEHe3 + fEHad', inplace=True)
    dataset.eval('fPLi = sqrt(fPxLi**2 + fPyLi**2 + fPzLi**2)', inplace=True)

    dataset.eval('fMassInvLi = sqrt(fELi**2 - fPLi**2)', inplace=True)
    dataset.query('fMassInvLi > 4.15314007', inplace=True)

    df_poshe3 = dataset.query('fSignHe3 > 0', inplace=False)
    df_neghe3 = dataset.query('fSignHe3 < 0', inplace=False)
    print(f'Positive He3: {df_poshe3.shape[0]}')
    print(f'Negative He3: {df_neghe3.shape[0]}')
