import numpy as np
import runDM

import matplotlib as mpl
import pylab as pl

font = { 'size'   : 16, 'family': 'serif'}
mpl.rcParams['xtick.major.size'] = 10
mpl.rcParams['xtick.major.width'] = 1
mpl.rcParams['xtick.minor.size'] = 4
mpl.rcParams['xtick.minor.width'] = 1
mpl.rcParams['ytick.major.size'] = 10
mpl.rcParams['ytick.major.width'] = 1
mpl.rcParams['ytick.minor.size'] = 4
mpl.rcParams['ytick.minor.width'] = 1
mpl.rc('font', **font)

#Z mass
mZ = 91.1875

#Nuclear energy scale
E_N = 1

c = sr.InitCouplings()
c[0] = 1.0
#c[-1] = 1.0

mV = np.logspace(0, 5, 100)

print "Check all this, because it doesn't make sense - check matrix ordering!"

pl.figure()

pl.plot(mV, np.dot(sr.EvolveMat(mV, E_N+ mV*0.0), c)[:,0], 'r-', linewidth=1.5, label=r'$c_u^{(V)}$')
pl.plot(mV, np.dot(sr.EvolveMat(mV, E_N+ mV*0.0), c)[:,1], 'b-',linewidth=1.5, label=r'$c_d^{(V)}$')
#pl.loglog(mV, np.dot(sr.EvolveMat(mV, E_N+ mV*0.0), c)[:,3], 'g-',linewidth=1.5, label=r'$c_s^{(V)}$')

pl.plot(mV, np.dot(sr.EvolveMat(mV, E_N+ mV*0.0), c)[:,8], 'r--',linewidth=1.5, label=r'$c_u^{(A)}$')
pl.plot(mV, np.dot(sr.EvolveMat(mV, E_N+ mV*0.0), c)[:,9], 'b--',linewidth=1.5, label=r'$c_d^{(A)}$')
pl.plot(mV, np.dot(sr.EvolveMat(mV, E_N+ mV*0.0), c)[:,11], 'g--', linewidth=1.5, label=r'$c_s^{(A)}$')

pl.legend(loc='lower left')

pl.axvline(mZ, linestyle='--', color='r')


pl.xscale('log')

pl.xlabel(r'$m_V$ [GeV]')
pl.ylabel(r'$c$')

pl.ylim(-2, 2)

pl.show()