# %load runDM.py
"""
Calculate low energy DM-SM couplings from high energy couplings.
"""

"""
Please contact Bradley Kavanagh (bradkav@gmail.com) for any questions,
problems, bugs and suggestions.
"""


import numpy as np
from scipy import linalg


#Define the Z-mass
mZ = 91.1875 #GeV

#Load in the evolution tables
gammaEMSM = np.loadtxt('data/gammaEMSM.dat')
U_match = np.loadtxt('data/Umatch.dat')
gammaSM = np.loadtxt('data/gammaSM.dat')


#%%
def initCouplings():
    """
    InitCouplings...
    """
    return np.zeros(16)


#%%
def setBenchmark(benchmark):
    """
    
    """
    if (benchmark == "Higgs"):
        return np.array([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0])
    elif (benchmark == "UniversalVector"):
        return np.array([1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.0])
    elif (benchmark == "UniversalAxial"):
        return np.array([-1.0,1.0,1.0,-1.0,1.0,-1.0,1.0,1.0,-1.0,1.0,-1.0,1.0,1.0,-1.0,1.0,0.0])
    elif (benchmark == "QuarksVector"):
        return np.array([1.0,1.0,1.0,0.0,0.0,1.0,1.0,1.0,0.0,0.0,1.0,1.0,1.0,0.0,0.0,0.0])
    elif (benchmark == "QuarksAxial"):
        return np.array([-1.0,1.0,1.0,0.0,0.0,-1.0,1.0,1.0,0.0,0.0,-1.0,1.0,1.0,0.0,0.0,0.0])
    elif (benchmark == "LeptonsVector"):
        return np.array([0.0,0.0,0.0,1.0,1.0,0.0,0.0,0.0,1.0,1.0,0.0,0.0,0.0,1.0,1.0,0.0])
    elif (benchmark == "LeptonsAxial"):
        return np.array([0.0,0.0,0.0,-1.0,1.0,0.0,0.0,0.0,-1.0,1.0,0.0,0.0,0.0,-1.0,1.0,0.0])
    elif (benchmark == "ThirdVector"):
        return np.array([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,1.0,1.0,1.0,1.0,0.0])
    elif (benchmark == "ThirdAxial"):
        return np.array([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,1.0,1.0,-1.0,1.0,0.0])
    
    else:
        print " Error in runDM.setBenchmark: Benchmark <<", benchmark, ">> not found..."
        print " Options are: 'Higgs', 'UniversalVector', 'UniversalAxial', 'QuarksVector', 'QuarksAxial', 'LeptonsVector', 'LeptonsAxial', 'ThirdVector', 'ThirdAxial'..."
        print " Returning empty coupling vector..."
        return np.zeros(16)
        

#Calculate evolution matrix from energy E1 to E2
#%%
def evolutionMat(E1, E2):
    """
    EvolveMat...
    """
    
    #Do bounds checking...
    
    #Check for E2 < 1
    #Precalculate the SMEM evolution matrix

    #Both energies below the Z-scale
    if ((E1 < mZ)and(E2 < mZ)):
        tnuclear = np.log(E2*1.0/E1)
        U_EMSM = linalg.expm(gammaEMSM*tnuclear)
        Emat = np.dot(U_EMSM,U_match)
    
    #Both energies above the Z-scale
    if ((E1 > mZ)and(E2 > mZ)):
        t = np.log(E1/E2)
        U_SM = linalg.expm(-gammaSM*t)
        Emat = np.dot(U_match,U_SM)
    
    #Energies across the Z-scale
    if ((E1 > mZ)and(E2 < mZ)):
        t = np.log(E1/mZ)
        tnuclear = np.log(E2/mZ)
        U_EMSM = linalg.expm(gammaEMSM*tnuclear)
        U_SM = linalg.expm(-gammaSM*t)
        Emat = np.dot(U_EMSM, np.dot(U_match, U_SM))
    
    return Emat
    
    
#%%
def evolveCouplings(c, E1, E2):
    """
    evolveCouplings...
    """
    
    #Need to check for size of c
    
    return np.dot(evolutionMat(E1, E2), c)


#%%
def lightqCouplings(c, E1, E2):
    """
    lightqCouplings...
    """
    #Might just want to fix E2...
    
    cf = evolveCouplings(c, E1, E2)
    return cf[[0,1,3,8,9,11]]
    
    