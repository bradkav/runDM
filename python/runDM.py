# %load runDM.py
"""
runDM version 1.0 [Python implementation]

Calculate low energy DM-SM couplings from high energy couplings,
taking into account RG evolution due to SM loops. See
arXiv:1605.04917 for further details.

Please contact Bradley Kavanagh (bradkav@gmail.com) for any questions,
problems, bugs and suggestions.
"""
from __future__ import print_function

import numpy as np
from scipy import linalg
from scipy.interpolate import interp1d
import sys

#-------------Initialisation---------------
#------------------------------------------

mZ = 91.1875 #Z-mass in GeV
mN = 0.938 #Nucleon mass in GeV

EvolutionSM = np.zeros([16,16,1401])
EvolutionEMSM = np.zeros([16,16,453])

#Pre-calculated values of t = Log[m_V/m_Z]
t_SM = np.linspace(0, 14, 1401)
t_EMSM = np.append(np.linspace(-4.51, 0.0,452)-0.003, 0)


#Load in the evolution tables
for i in range(1401):
    #Truncate to a few decimal places, to prevent rounding errors in the filenames
    s = str(np.around(t_SM[i], 3))
    if (t_SM[i] == int(t_SM[i])):
        s = str(int(t_SM[i]))
    EvolutionSM[:,:,i] = np.loadtxt('data/EvolutionSM_t=' + s + '.dat')

for i in range(453):
    #Truncate to a few decimal places, to prevent rounding errors in the filenames
    s = str(np.around(t_EMSM[i], 3))
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

#------------------------------------------


#%% Initialise empty coupling vector
def initCouplings():
    """
    initCouplings()
    Returns a numpy array with 16 elements, all set to zero,
    for use in initialising coupling vectors.
    """
    return np.zeros(16)


#%% Generate coupling vectors with preset operator structures
def setBenchmark(benchmarkID):
    """
    setBenchmark(benchmarkID)
    Returns a numpy array with 16 elements, corresponding to the
    vector of couplings defined in Eq. 4 of the runDM manual. 
    The value of the couplings is defined by the string benchmarkID.
    
    Possible choices for benchmarkID are:
        'Higgs', 'UniversalVector', 'UniversalAxial',
        'QuarksVector', 'QuarksAxial', 'LeptonsVector',
        'LeptonsAxial', 'ThirdVector', 'ThirdAxial'
    """
    if (benchmarkID == "Higgs"):
        return np.array([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0])
    elif (benchmarkID == "UniversalVector"):
        return np.array([1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.0])
    elif (benchmarkID == "UniversalAxial"):
        return np.array([-1.0,1.0,1.0,-1.0,1.0,-1.0,1.0,1.0,-1.0,1.0,-1.0,1.0,1.0,-1.0,1.0,0.0])
    elif (benchmarkID == "QuarksVector"):
        return np.array([1.0,1.0,1.0,0.0,0.0,1.0,1.0,1.0,0.0,0.0,1.0,1.0,1.0,0.0,0.0,0.0])
    elif (benchmarkID == "QuarksAxial"):
        return np.array([-1.0,1.0,1.0,0.0,0.0,-1.0,1.0,1.0,0.0,0.0,-1.0,1.0,1.0,0.0,0.0,0.0])
    elif (benchmarkID == "LeptonsVector"):
        return np.array([0.0,0.0,0.0,1.0,1.0,0.0,0.0,0.0,1.0,1.0,0.0,0.0,0.0,1.0,1.0,0.0])
    elif (benchmarkID == "LeptonsAxial"):
        return np.array([0.0,0.0,0.0,-1.0,1.0,0.0,0.0,0.0,-1.0,1.0,0.0,0.0,0.0,-1.0,1.0,0.0])
    elif (benchmarkID == "ThirdVector"):
        return np.array([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,1.0,1.0,1.0,1.0,0.0])
    elif (benchmarkID == "ThirdAxial"):
        return np.array([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,1.0,1.0,-1.0,1.0,0.0])
    
    else:
        print(" Error in runDM.setbenchmark: benchmarkID <<", benchmarkID, ">> not found...")
        print(" Options are: 'Higgs', 'UniversalVector', 'UniversalAxial', 'QuarksVector', 'QuarksAxial', 'LeptonsVector', 'LeptonsAxial', 'ThirdVector', 'ThirdAxial'...")
        print(" Returning empty coupling vector...")
        return np.zeros(16)
        
    
    
#%% Calculate matrix to evolve couplings
def evolutionMat(E1, E2):
    """
    evolutionMat(E1, E2)
    Calculate 16x16 matrix for evolving coupling vector
    between energies E1 and E2 (in GeV). evolutionMat
    takes care of the relative values of E1 and E2.
    
    Requires E1, E2 in range [1, 1e8] GeV.
    
    Note that evolutionMat is NOT vectorized - it can only
    accept floats for E1 and E2.
    
    Input:
        E1 - energy to start from (in GeV)
        E2 - energy to run to (in GeV)
    
    Output:
        Returns a 16x16 numpy array containing the evolution matrix
    """
    
    #Check to see if E1 or E2 is a list
    if ((hasattr(E1, "__len__"))or(hasattr(E2,"__len__"))):
        sys.exit(" Error in runDM.evolutionMat: E1 and E2 must both be floats")
                
    
    t1 = np.log(E1/mZ)
    t2 = np.log(E2/mZ)
    
    #Check ranges of t1 and t2
    if not(np.log(1.0/mZ) <= t1 <= 14.0):
        sys.exit(" Error in runDM.evolutionMat: E1 out of range. Require 1 GeV <= E1 <= 1e8 GeV")

    if not(np.log(1.0/mZ) <= t2 <= 14.0):
        sys.exit(" Error in runDM.evolutionMat: E2 out of range. Require 1 GeV <= E2 <= 1e8 GeV")

    #Both energies below the Z-mass
    if ((t1 <= 0)and(t2 <= 0)):
        Emat = np.dot(np.dot(linalg.inv(UevolutionBELOW(t2)), UevolutionBELOW(t1)),Umatch)
    
    #Both energies above the Z-mass
    if ((t1 >= 0)and(t2 >= 0)):
        Emat = np.dot(linalg.inv(UevolutionABOVE(t2)), UevolutionABOVE(t1))
    
    #Energies either side of the Z-mass
    if ((t1 >= 0)and(t2 <= 0)):
        EmatBELOW = np.dot(linalg.inv(UevolutionBELOW(t2)), UevolutionBELOW(0))
        Emat = np.dot(EmatBELOW, np.dot(Umatch, UevolutionABOVE(t1))) 
    
    #Energies either side of Z-mass (in wrong order)
    if ((t1 < 0)and(t2 > 0)):
        sys.exit(" Error in runDM.evolutionMat: E1 < mZ, E2 > mZ not supported - matching is not unique.")
    
    return Emat    
    

    
#%% Evolve couplings between E1 and E2
def runCouplings(c, E1, E2):
    """
    runCouplings(c, E1, E2)
    Calculate running of couplings c between two energies
    E1 and E2. If E2 > mZ, the output is an array of 
    couplings in the EW-unbroken phase (Eq. 4 of the manual). 
    If E2 < mZ, the output is an array of couplings in
    the EW-broken phase (Eq. 6 of the manual).
    
    Note that E1 < mZ with E2 > mZ is not allowed.
    
    Input:
        c - numpy array with 16 elements, with values corresponding
            to those defined in Eq. 4 of the runDM manual
        E1 - energy (in GeV) at which c is defined. E1 may be
            a scalar or a 1-d numpy array.
        E2 - energy (in GeV) to run to. E2 may be a scalar or 1-d
            numpy array (with same length as E1).
    
    Output:
        Returns array with length 16 in the final dimension
        (corresponding either to Eq. 4 or Eq. 6 of the manual).
        The full dimensions of the array will be (Len(E1), 16).
    """
    
    #Check length of coupling vector c
    if not(hasattr(c, "__len__")):
        sys.exit(" Error in runDM.runCouplings: c must be an array with 16 elements")
    if (len(c) != 16):
        sys.exit(" Error in runDM.runCouplings: c must be an array with 16 elements")
    
    
    #If E1 and E2 are scalar
    if not((hasattr(E1, "__len__"))or(hasattr(E2,"__len__"))):
        return np.dot(evolutionMat(E1, E2), c)
       
    #If E1 or E2 are arrays, need to check to make sure correct
    #array dimensions are returned 
    
    #Both E1, E2 are arrays (check they are same length...)
    if ((hasattr(E1, "__len__"))and(hasattr(E2, "__len__"))):
        n1 = len(E1)
        if (len(E2) != n1):
            sys.exit(" Error in runDM.runCouplings: E1 and E2 must have same length (or be scalar)")
        else:
            result = np.zeros([n1,16])
            for i in range(n1):
                result[:, i] = np.dot(evolutionMat(E1[i], E2[i]), c)

    #Only E1 is an array
    elif (hasattr(E1, "__len__")):
        n1 = len(E1)
        result = np.zeros([n1,16])
        for i in range(n1):
            result[:, i] = np.dot(evolutionMat(E1[i], E2), c)

    #Only E2 is an array
    elif (hasattr(E2, "__len__")):
        n2 = len(E2)
        result = np.zeros([n2,16])
        for i in range(n2):
            result[i, :] = np.dot(evolutionMat(E1, E2[i]), c)


    return result



#%% Calculate couplings to light quarks at the nuclear scale
def DDCouplingsQuarks(c, E1):
    """
    DDCouplingsQuarks(c, E1)
    Calculate vector (V) and axial-vector (A) couplings
    to u, d, s quarks at the nuclear energy scale
    starting from the high energy coupling vector c,
    defined at energy E1 (in GeV).

    Input:
        c - numpy array with 16 elements, with values corresponding
            to those defined in Eq. 4 of the runDM manual
        E1 - energy (in GeV) at which c is defined. E1 may be
            a scalar or a 1-d numpy array.
    
    Output:
        Returns array with length 5 in the final dimension,
        corresponding to (CVu, CVd, CAu, CAd, CAs).
        The dimensions of the array will be (Len(E1), 5).
    """
    
    #Check length of coupling vector c
    if not(hasattr(c, "__len__")):
        sys.exit(" Error in runDM.runCouplings: c must be an array with 16 elements")
    if (len(c) != 16):
        sys.exit(" Error in runDM.runCouplings: c must be an array with 16 elements")
    
    #If E1 is scalar
    if not(hasattr(E1, "__len__")):
        return runCouplings(c, E1, 1.0)[[0,1,8,9,11]]

    #If E1 is an array
    else:
        n1 = len(E1)
        result = np.zeros([n1,5])
        for i in range(n1):
            result[i,:] = runCouplings(c, E1[i], 1.0)[[0,1,8,9,11]]
            
        return result
    
#%% Calculate non-relativistic couplings to protons
def DDCouplingsProton(c, E1, mx, DMcurrent):
    
	#From quarks to nucleons
	#Values from arXiv:1202.1292
    deltau_p = 0.84
    deltad_p = -0.44
    deltas_p = -0.03
    
    #Get quark couplings
    cuV, cdV, cuA, cdA, csA = DDCouplingsQuarks(c, E1)

    #Calculate non-relativistic proton couplings
    #Note the number of operators is shifted by 1
    #because python uses zero-indexed arrays
    lambda_p = np.zeros(12)
    if (DMcurrent == "scalar"):
        lambda_p[0] = 4*mx*mN*(2*cuV+ cdV)
        lambda_p[6] = -8*mx*mN*(deltau_p*cuA + deltad_p*cdA + deltas_p*csA)
    elif (DMcurrent == "vector"):
        lambda_p[0] = 4*mx*mN *(2*cuV + cdV)
        lambda_p[6] = -8*mx*mN*(deltau_p*cuA + deltad_p*cdA + deltas_p*csA)
        lambda_p[8] = 8*mN*(deltau_p*cuA + deltad_p*cdA + deltas_p*csA)
    elif (DMcurrent == "axial-vector"):
        lambda_p[3] = -16*mx*mN*(deltau_p*cuA + deltad_p*cdA + deltas_p*csA)
        lambda_p[7] = 8*mx*mN*(2*cuV + cdV)
        lambda_p[8] = 8*mx*(2*cuV + cdV)

    return lambda_p/E1**2
    
#%% Calculate non-relativistic couplings to neutrons
def DDCouplingsNeutron(c, E1, mx, DMcurrent):

    #From quarks to nucleons
    #Values from arXiv:1202.1292
    deltau_p = 0.84
    deltad_p = -0.44
    deltas_p = -0.03

    #Get quark couplings
    cuV, cdV, cuA, cdA, csA = DDCouplingsQuarks(c, E1)

    #Calculate non-relativistic neutron couplings
    #Note the number of operators is shifted by 1
    #because python uses zero-indexed arrays
    lambda_n = np.zeros(12)
    if (DMcurrent == "scalar"):
        lambda_n[0] = 4*mx*mN*(cuV + 2*cdV)
        lambda_n[6] = -8*mx*mN*(deltad_p*cuA + deltau_p*cdA + deltas_p*csA)
    elif (DMcurrent == "vector"):
        lambda_n[0] = 4*mx*mN*(cuV + 2*cdV)
        lambda_n[6] = -8*mx*mN*(deltad_p*cuA + deltau_p*cdA + deltas_p*csA)
        lambda_n[8] = 8*mN*(deltad_p*cuA + deltau_p*cdA + deltas_p*csA)
    elif (DMcurrent == "axial-vector"):
        lambda_n[3] = -16*mx*mN*(deltad_p*cuA + deltau_p*cdA + deltas_p*csA)
        lambda_n[7] = 8*mx*mN*(cuV + 2*cdV)
        lambda_n[8] = 8*mx*(cuV + 2.0*cdV)
    
    return lambda_n/E1**2


#%% Calculate non-relativistic couplings to nucleons
def DDCouplingsNR(c, E1, mx, DMcurrent, N):
    """
    DDCouplingsNR(c, E1, mx, DMcurrent, N)
	Calculate coefficients of the non-relativistic DM-nucleon
	operators at the nuclear energy scale, with numbering as 
	in arXiv:1307.5955, taking into account running of the 
	couplings from high energy E1 (in GeV). The result will 
	depend on the structure of the DM interaction and on the 
    DM mass. 

    Input:
        c - vector with 16 elements, with values corresponding
            to those defined in Eq. 4 of the runDM manual
        E1 - energy (in GeV) at which c is defined. E1 may be
            a scalar or an array.
		mx - DM mass in GeV. mx must be a single number.
		DMcurrent - string specifying the DM interaction current.
			The options are 'scalar', 'vector' and 'axial-vector'.
			The corresponding definitions of the DM currents are
			given in Eq. 8 of the manual.
		N - string specifying the nucleon type: 'p' for proton,
			'n' for neutron.

    Output:
        Returns array with length 12 in the final dimension,
        corresponding to the coefficients of the first 12
		non-relativistic DM-nucleon operators listed in 
		arXiv:1307.5955. The coefficients include a factor
        of 1/E1^2.  The dimensions of the array will
		be (Len(E1), 12)."
    """
    
    #Check parameter values
    if N not in ("p", "n"):
        sys.exit(" Error in runDM.DDCouplingsNR: N must be either 'p' or 'n'.")
    if DMcurrent not in ("scalar", "vector", "axial-vector"):
        sys.exit(" Error in runDM.DDCouplingsNR: DMcurrent must be 'scalar', 'vector' or 'axial-vector'.")
    
    #Check length of coupling vector c
    if not(hasattr(c, "__len__")):
        sys.exit(" Error in runDM.DDCouplingsNR: c must be an array with 16 elements")
    if (len(c) != 16):
        sys.exit(" Error in runDM.DDCouplingsNR: c must be an array with 16 elements")
    
    #Check that mx is a single number
    if hasattr(mx, "__len__"):
        sys.exit(" Error in runDM.DDCouplingsNR: mx must be a single number.")

    #Determine nucleon type
    if (N == "p"):
        DDfunc = lambda E: DDCouplingsProton(c, E, mx, DMcurrent)
    elif (N == "n"):
        DDfunc = lambda E: DDCouplingsNeutron(c, E, mx, DMcurrent)
    
    #If E1 is scalar
    if not(hasattr(E1, "__len__")):
        return DDfunc(E1)

    #If E1 is an array
    else:
        n1 = len(E1)
        result = np.zeros([n1,12])
        for i in range(n1):
            result[i,:] = DDfunc(E1[i])
        return result
        