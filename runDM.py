# %load runDM.py
"""
Calculate low energy DM-SM couplings from high energy couplings.
"""


import numpy as np
from scipy import linalg


#Define the Z-mass
mZ = 91.1875 #GeV

#Load in the evolution tables
gammaEMSM = np.loadtxt('gammaEMSM.dat')
U_match = np.loadtxt('Umatch.dat')
gammaSM = np.loadtxt('gammaSM.dat')


#Calculate evolution matrix from energy E1 to E2
#%%
def EvolveMat(E1list, E2list):
    """
    EvolveMat...
    """

    Emat = np.zeros([len(E1list),16,16])
    
    for i in range(len(E1list)):
    
        E1 = E1list[i]
        E2 = E2list[i]
        
        #Both energies below the Z-scale
        if ((E1 < mZ)and(E2 < mZ)):
            tnuclear = np.log(E2*1.0/E1)
            U_EMSM = linalg.expm(gammaEMSM*tnuclear)
            Emat[i,:,:] = np.dot(U_EMSM,U_match)
        
        #Both energies above the Z-scale
        if ((E1 > mZ)and(E2 > mZ)):
            t = np.log(E1/E2)
            U_SM = linalg.expm(-gammaSM*t)
            Emat[i,:,:] = np.dot(U_match,U_SM)
        
        #Energies across the Z-scale
        if ((E1 > mZ)and(E2 < mZ)):
            t = np.log(E1/mZ)
            tnuclear = np.log(E2/mZ)
            U_EMSM = linalg.expm(gammaEMSM*tnuclear)
            U_SM = linalg.expm(-gammaSM*t)
            Emat[i,:,:] = np.dot(U_EMSM, np.dot(U_match, U_SM))
    
    return Emat

#%%
def InitCouplings():
    """
    InitCouplings...
    """
    return np.zeros(16)


