'''
    Coulomb CF implementation
'''

import numpy as np
from functools import partial
from ROOT import TF1, TFile

M_PI = 139.57039 # MeV
M_PR = 938.272081 # MeV
M_HE3 = 2808.391 # MeV
MU_PRHE = 1 / (1/M_PR + 1/M_HE3) # MeV
HBAR_C = 197.327 # MeV fm
ALPHA_EM = 1/137.036

def KGamow(x: np.ndarray, pars: np.ndarray) -> float:
    '''
        Coulomb correlation function for two charged pions using a point like source
        q = x[0]: relative momentum in MeV
        z1 = pars[0]: charge of particle 1
        z2 = pars[1]: charge of particle 2
        mu = pars[2]: reduced mass in MeV
    '''

    z1 = pars[0]
    z2 = pars[1]
    mu = pars[2]
    eta = ALPHA_EM * mu * z1 * z2 / (x[0])
    return 2*np.pi*eta/(np.exp(2*np.pi*eta) - 1)


'''
    Coulomb correction for a correlation function using an extended Levy source
    Gaussian case is a special case of the Levy source with alpha = 2
'''
aA, aB, aC, aD, aE = 0.36060, -0.54508, 0.03475, -1.30389, 0.00378
bA, bB, bC, bD, bE, bF = 2.04017, 0.55972, 2.47224, -1.26815, -0.11767, 0.52738
cA, cB, cC, cD, cE, cF = -1.00015, 0.00012, 0.00008, 0.26986, 0.00003, 1.75202
dA, dB, dC, dD, dE, dF = 0.00263, -0.13124, -0.83149, 1.57528, 0.27568, 0.04937

def A(alpha, R, aA, aB, aC, aD, aE):
    return (aA * alpha + aB) ** 2 + (aC * R + aD) ** 2 + aE * (alpha * R + 1) ** 2
def B(alpha, R, bA, bB, bC, bD, bE, bF):
    return 1 + bA * R ** bB - alpha * bC / (alpha ** 2 * R * (alpha * bD + bE * R ** bF))
def C(alpha, R, cA, cB, cC, cD, cE, cF):
    return (cA + alpha * cB + cC * R ** cD) / cE * (alpha / R) ** cF
def D(alpha, R, dA, dB, dC, dD, dE, dF):
    return dA + (R ** dB + dC * alpha ** dF)/(R ** dD * alpha ** dE)

def KmodGeneral(q, alpha, R, z1, z2, mu, aA, aB, aC, aD, aE, bA, bB, bC, bD, bE, bF, cA, cB, cC, cD, cE, cF, dA, dB, dC, dD, dE, dF):
    return (1 + A(alpha, R, aA, aB, aC, aD, aE) * 
            ALPHA_EM * z1 * z2 * mu * np.pi * R / (alpha * HBAR_C) / (1 +
            B(alpha, R, bA, bB, bC, bD, bE, bF) * (q * R / (alpha * HBAR_C)) +
            C(alpha, R, cA, cB, cC, cD, cE, cF) * (q * R / (alpha * HBAR_C)) ** 2 +
            D(alpha, R, dA, dB, dC, dD, dE, dF) * (q * R / (alpha * HBAR_C)) ** 4))

KmodGeneralFixed = partial(KmodGeneral, aA=aA, aB=aB, aC=aC, aD=aD, aE=aE, bA=bA, bB=bB, bC=bC, bD=bD, bE=bE, bF=bF, cA=cA, cB=cB, cC=cC, cD=cD, cE=cE, cF=cF, dA=dA, dB=dB, dC=dC, dD=dD, dE=dE, dF=dF)

def KmodCMS(q, R, z1, z2, mu):
    return (1 + (ALPHA_EM * np.pi * z1 * z2 * mu * R) / (1.26 * HBAR_C + q * R))

def KLevyGeneral(x: np.ndarray, pars: np.ndarray) -> float:
    '''
        Coulomb correction for a correlation function using an extended Levy source
        Gaussian case is a special case of the Levy source with alpha = 2
        k = x[0]: k* in MeV
        alpha = pars[0]: Levy source parameter
        R = pars[1]: Source size in fm
        z1 = pars[2]: charge of particle 1
        z2 = pars[3]: charge of particle 2
        mu = pars[4]: reduced mass in MeV
    '''
    q = 2*x[0]
    alpha = pars[0]
    R = pars[1]
    z1 = pars[2]
    z2 = pars[3]
    mu = pars[4]
    aA, aB, aC, aD, aE = pars[5], pars[6], pars[7], pars[8], pars[9]
    bA, bB, bC, bD, bE, bF = pars[10], pars[11], pars[12], pars[13], pars[14], pars[15]
    cA, cB, cC, cD, cE, cF = pars[16], pars[17], pars[18], pars[19], pars[20], pars[21]
    dA, dB, dC, dD, dE, dF = pars[22], pars[23], pars[24], pars[25], pars[26], pars[27]
    return KGamow([q], [z1, z2, mu]) * KmodGeneral(q, alpha, R, z1, z2, mu, aA, aB, aC, aD, aE, bA, bB, bC, bD, bE, bF, cA, cB, cC, cD, cE, cF, dA, dB, dC, dD, dE, dF)

def KLevyGeneralFixed(x: np.ndarray, pars: np.ndarray) -> float:
    '''
        Coulomb correction for a correlation function using an extended Levy source
        Gaussian case is a special case of the Levy source with alpha = 2
        k = x[0]: k* in MeV
        alpha = pars[0]: Levy source parameter
        R = pars[1]: Source size in fm
        z1 = pars[2]: charge of particle 1
        z2 = pars[3]: charge of particle 2
        mu = pars[4]: reduced mass in MeV
    '''
    q = 2*x[0]
    alpha = pars[0]
    R = pars[1]
    z1 = pars[2]
    z2 = pars[3]
    mu = pars[4]
    return KGamow([q], [z1, z2, mu]) * KmodGeneralFixed(q, alpha, R, z1, z2, mu)

def KLevyCMS(x: np.ndarray, pars: np.ndarray) -> float:
    '''
        Coulomb correction for a correlation function using an extended Levy source
        Gaussian case is a special case of the Levy source with alpha = 2
        k = x[0]: k* in MeV
        alpha = pars[0]: NOT USED
        R = pars[1]: Source size in fm
        z1 = pars[2]: charge of particle 1
        z2 = pars[3]: charge of particle 2
        mu = pars[4]: reduced mass in MeV
    '''
    q = 2*x[0]
    R = pars[1]
    z1 = pars[2]
    z2 = pars[3]
    mu = pars[4]
    return KGamow([q], [z1, z2, mu]) * KmodCMS(q, R, z1, z2, mu)


def fit_cf(infile_path, klevy, outfile, name):
    infile = TFile.Open(infile_path)
    cf = infile.Get("hHe3_p_Coul_CF_LS")
    #for ibin in range(1, cf.GetNbinsX() + 1):
    #    cf.SetBinError(ibin, 0.05)
    cf.GetXaxis().SetLimits(cf.GetXaxis().GetXmin() * 1000, cf.GetXaxis().GetXmax() * 1000)
    cf.Fit(klevy, "RMS+")
    print(f"{name}: {klevy.GetChisquare()}/{klevy.GetNDF()}")
    outfile.cd()
    cf.Write(name)

if __name__ == "__main__":
    outfile = TFile.Open("/home/galucia/antiLithium4/femto/output/coulomb.root", "RECREATE")
    
    '''
    klevy = TF1("klevy", KLevyGeneral, 2, 400, ndim=1, npar=28)
    for ipar, par in zip(range(5, 28), [aA, aB, aC, aD, aE, bA, bB, bC, bD, bE, bF, cA, cB, cC, cD, cE, cF, dA, dB, dC, dD, dE, dF]):
        klevy.SetParameter(ipar, par)
    klevy.FixParameter(0, 2)
    klevy.FixParameter(1, 6)
    klevy.FixParameter(2, 1) # z1 = 1; p
    klevy.FixParameter(3, 2) # z2 = 2; He3
    klevy.FixParameter(4, MU_PRHE) # mu = reduced mass of p and He3
    klevy.Write()

    infile_path0_10 = '/home/galucia/antiLithium4/analysis/output/CATS/CATS_cent0_10.root'
    fit_cf(infile_path0_10, klevy, outfile, '0_10')
    '''

    klevy = TF1("klevy", KLevyCMS, 2, 400, ndim=1, npar=5)
    klevy.FixParameter(0, 2)
    klevy.FixParameter(1, 6)
    klevy.FixParameter(2, 1) # z1 = 1; p
    klevy.FixParameter(3, 2) # z2 = 2; He3
    klevy.FixParameter(4, MU_PRHE) # mu = reduced mass of p and He3
    klevy.Write()

    infile_path0_10 = '/home/galucia/antiLithium4/analysis/output/CATS/CATS_cent0_10.root'
    fit_cf(infile_path0_10, klevy, outfile, '0_10')

    '''
    klevy = TF1("klevy", KLevyGeneralFixed, 2, 200, ndim=1, npar=5)
    klevy.FixParameter(0, 2)
    klevy.SetParameter(1, 4)
    klevy.FixParameter(2, 1) # z1 = 1; p
    klevy.FixParameter(3, 2) # z2 = 2; He3
    klevy.FixParameter(4, MU_PRHE) # mu = m_pi
    klevy.Write()
    kgamow = TF1("kgamow", KGamow, 2, 800, ndim=1, npar=3)
    kgamow.SetParameters(1, 2, MU_PRHE)
    kgamow.Write()

    klevy.SetParLimits(1, 0, 10)
    infile_path0_10 = '/home/galucia/antiLithium4/analysis/output/CATS/CATS_cent0_10.root'
    klevy.SetParameter(1, 6)
    fit_cf(infile_path0_10, klevy, outfile, '0_10')

    infile_path10_30 = '/home/galucia/antiLithium4/analysis/output/CATS/CATS_cent10_30.root'
    klevy.SetParameter(1, 4.5)
    fit_cf(infile_path10_30, klevy, outfile, '10_30')

    infile_path30_50 = '/home/galucia/antiLithium4/analysis/output/CATS/CATS_cent30_50.root'
    klevy.SetParameter(1, 3)
    fit_cf(infile_path30_50, klevy, outfile, '30_50')
    '''

    outfile.Close()
