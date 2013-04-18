import pylab as p
import numpy as np
from numpy import exp, cosh
from scipy.optimize import brentq, leastsq
from parse_overlap import OverlapPoint, OverlapWall

def JKsigma(cfnc):
    '''Errorbars given correlator average and jackknife correlators.'''
    ave, JKvals = cfnc[0], cfnc[1:]
    N = len(JKvals)
    diffs = (ave-JKvals)*(ave-JKvals)
    return np.sqrt(np.sum(diffs, axis=0)*(1-1./N))

def plot_correlator(cfnc):
    ave, sigma = cfnc[0], JKsigma(cfnc)
    p.figure()
    p.xlabel('$t$')
    p.ylabel('$C[t]$')
    p.yscale('log')
    p.errorbar(range(96), ave.real,sigma.real, fmt='k-')
    p.show()
    
def fit_twopoint(xarr, yarr, earr, T=96.): #might add t1, t2 here
    '''Fit two-point correlator to the cosh form.'''
    fitfunc = lambda p, x: p[0]*(2*exp(-p[1]*T/2.)*cosh(p[1]*(x-T/2.)))
    errfunc = lambda p, x, y, err: (y - fitfunc(p,x))/err
    
    p0 = np.array([0.5, 0.5])  # Initial guess.
    p1, success = leastsq(errfunc, p0, args=(xarr,yarr,earr),
                                   full_output=0)
    return p1

    
def fold(cfnc):
    '''Fold correlators around the mid-point.'''
    T = len(cfnc[0])
    assert (T % 2 == 0)  
    frow, c = cfnc[:,:1], cfnc[:,1:]  # Split off separation=0 point. 
    c = (c + np.fliplr(c))/2.  # Mid-point is fixed.
    return np.hstack((frow, c))  # (Second half is now redundant.)
    
def m_eff(cfnc, T):
    '''Fit m_eff using the cosh form. 
    
    Pass in a -- single -- correlation function.
    T is passed explicitly since we may not be dealing with full correlator.
    '''
    f = lambda m, t, A: A - np.cosh(m*(t-T/2.))/np.cosh(m*(t+1-T/2.))
    return [brentq(f, 0, 10, args=(t,cfnc[t]/cfnc[t+1]))
            for t in range(1, T/2-1)]  # Works better than fsolve.
    
def cosh_effmass(cfnc, foldQ=True):
    '''Compute effective mass of correlation functions, using cosh form.'''
    T = len(cfnc[0])
    if foldQ:
        cfnc = fold(cfnc.real)
    return np.array([m_eff(c, T) for c in cfnc])

    
def naive_effmass(cfnc, foldQ=False):
    '''Compute naive effective mass of correlation functions.'''
    if foldQ:
        cfnc = fold(cfnc)
    cplus = np.roll(cfnc, 1, axis=1).real
    return np.abs(np.log(cfnc/cplus))
    
def plot_effmass2(cfnc):
    effmass = cosh_effmass(cfnc)
    ave, sigma = effmass[0], JKsigma(effmass)
    p.figure()
    p.xlabel('$t$')
    p.ylabel('$|\\frac{C[t+1]}{C[t]}|$')
    p.errorbar(range(len(ave)), ave, sigma, fmt='ko') # Fix.
    p.show()
    
def plot_effmass(cfnc):
    effmass = naive_effmass(cfnc)
    ave, sigma = effmass[0], JKsigma(effmass)
    p.figure()
    p.xlabel('$t$')
    p.ylabel('$|\\frac{C[t+1]}{C[t]}|$')
    p.errorbar(range(96), ave, sigma, fmt='ko')
    p.show()
    
def main():
    ms = 0.0495
    wall = OverlapWall(0.55, ms)
    point = OverlapPoint(0.551, ms)
    
    plot_correlator(wall.pscalar)
#    plot_correlator(point.pscalar)
#    plot_correlator(wall.vector)
#    plot_correlator(point.vector)
#    plot_effmass(wall.pscalar)
#    plot_effmass2(wall.pscalar)
#    plot_effmass(point.pscalar)
#    plot_effmass(wall.vector)
#    plot_effmass(point.vector)
    
if __name__ == "__main__":
    main()
