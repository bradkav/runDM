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
from scipy.interpolate import interp1d


#Define the Z-mass
mZ = 91.1875 #GeV

EvolutionSM = np.zeros([16,16,1401])
EvolutionEMSM = np.zeros([16,16,453])

t_SM = np.linspace(0, 14, 1401)
t_EMSM = np.append(np.linspace(-4.51, 0.0,452)-0.003, 0)


#Load in the evolution tables
for i in range(1401):
    s = str(t_SM[i])
    if (t_SM[i] == int(t_SM[i])):
        s = str(int(t_SM[i]))
    EvolutionSM[:,:,i] = np.loadtxt('data/EvolutionSM_t=' + s + '.dat')

for i in range(453):
    s = str(t_EMSM[i])
    if (t_EMSM[i] == int(t_EMSM[i])):
        s = str(int(t_EMSM[i]))
    EvolutionEMSM[:,:,i] = np.loadtxt('data/EvolutionEMSM_t=' + s + '.dat')

Umatch = np.loadtxt('data/Umatch.dat')

#Correct the value of t_EMSM slightly
# (because Log(1/mZ)~-4.51292)
t_EMSM[:-1] = t_EMSM[:-1]+0.00008


#Define interpolating functions
UevolutionABOVE = interp1d(t_SM, EvolutionSM)
UevolutionBELOW = interp1d(t_EMSM, EvolutionEMSM) 


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
def evolutionMat0(E1, E2):
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
    
#Calculate evolution matrix from energy E1 to E2
#%%
def evolutionMat(E1, E2):
    """
    EvolveMat...
    """
    
    t1 = np.clip(np.log(E1/mZ), np.log(1.0/mZ), 14.0)
    t2 = np.clip(np.log(E2/mZ), np.log(1.0/mZ), 14.0)
    
    #Do bounds checking...
    
    #Check for E2 < 1
    #Precalculate the SMEM evolution matrix

    #Both energies below the Z-scale
    if ((E1 < mZ)and(E2 < mZ)):
        Emat = np.dot(np.dot(linalg.inv(UevolutionBELOW(t2)), UevolutionBELOW(t1)),Umatch)
    
    #Both energies above the Z-scale
    if ((E1 > mZ)and(E2 > mZ)):
        Emat = np.dot(linalg.inv(UevolutionABOVE(t2)), UevolutionABOVE(t1))
    
    #Energies across the Z-scale
    if ((E1 > mZ)and(E2 < mZ)):
        EmatBELOW = np.dot(linalg.inv(UevolutionBELOW(t2)), UevolutionBELOW(0))
        Emat = np.dot(EmatBELOW, np.dot(Umatch, UevolutionABOVE(t1))) 
    
    return Emat    
    
#%%
def evolveCouplings(c, E1, E2):
    """
    evolveCouplings...
    """
    
    #Need to check for size of c
    
    return np.dot(evolutionMat(E1, E2), c)


#%%
def DDCouplings(c, E1):
    """
    lightqCouplings...
    """
    #Might just want to fix E2...
    
    cf = evolveCouplings(c, E1, 1.0)
    return cf[[0,1,8,9,11]]
    
    