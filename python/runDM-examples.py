"""
Example script for the python implementation of runDM. See also the 
runDM-examples.ipynb ipython notebook for more detailed explanation.

Please contact Bradley Kavanagh (bradkav@gmail.com) for any questions,
problems, bugs and suggestions.
"""

import numpy as np
import matplotlib
from matplotlib import pyplot as pl

#Import the runDM module
import runDM

#First, let's specify the high-energy couplings. This will be an 1-D array with 16 elements. runDM comes with a number of pre-defined benchmarks, which can be accessed using setBenchmark.

c_high = runDM.setBenchmark("UniversalVector")
print "Vector coupling to all SM fermions:", c_high, "\n"

c_high = runDM.setBenchmark("QuarksAxial")
print "Axial-vector coupling to all quarks:", c_high, "\n"


#Alternatively, you can specify each coupling individually. You can use InitCouplings() to generate an empty array of couplings and then go ahead. But any array of 16 elements with do.

c_high = runDM.initCouplings()
c_high[0] = 1.0
c_high[1] = -1.0
c_high[12] = 1.0
print "User-defined couplings:", c_high, "\n"

#From these high energy couplings (defined at some energy E_1), you can obtain the couplings at a different energy scale E_2 by using runCouplings(c, E_1, E_2). See the manual for more details on runCouplings.

E1 = 1000
E2 = 1
c_low = runDM.runCouplings(c_high, E1, E2)
print "Low energy couplings:", c_low, "\n"

#If we're only interested in direct detection experiments, we can use the function DDCouplings(c, E_1). In this case, the code evolves the couplings from energy $E_1$, down to the nuclear energy scale ~ 1 GeV. The output is an array with 5 elements, the vector and axial-vector couplings to the light quarks.

c_q = runDM.DDCouplings(c_high, E1)

couplings_str = ['c_V^u','c_V^d','c_A^u','c_A^d','c_A^s']

for k in range(5):
    print couplings_str[k], "=", c_q[k]
print " "

#Finally, let's take a look at the value of the low-energy light quark couplings (evaluated at mu_N ~ 1 GeV) as a function of the mediator mass m_V. 


#Set the value of the high energy couplings
c_high = runDM.setBenchmark("QuarksAxial")

#Calculate the low energy couplings
mV = np.logspace(0, 6, 1000)
c_q = runDM.DDCouplings(c_high, mV)

#Now let's do some plotting
f, axarr = pl.subplots(3,2 ,figsize=(8,8))

for k in range(5):
    if (k < 2): #Vector currents
        ax = axarr[k%3, 0]
    else:       #Axial-vector currents
        ax = axarr[(k+1)%3, 1]
        
    ax.semilogx(mV, c_q[:,k])
    ax.set_xlabel(r'$m_V$ [GeV]', fontsize=18.0)
    ax.set_ylabel(r'$'+couplings_str[k]+'$', fontsize=20.0)
    ax.axvline(91.1875, color='k', linestyle='--')
    ax.set_xlim(1.0, 10**6)
    ax.get_yticklabels()[-1].set_visible(False)
    ax.tick_params(axis='both', labelsize=12.0)
    
axarr[2,0].set_axis_off()
pl.tight_layout()
pl.show()