import pylab as p
import numpy as np
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
    
def naive_effmass(cfnc):
    '''Compute effective mass of correlation function.'''
    cplus = np.roll(cfnc, 1, axis=1).real
    return np.abs(np.log(cfnc/cplus))
    
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
    plot_correlator(point.pscalar)
    plot_correlator(wall.vector)
    plot_correlator(point.vector)
#    plot_effmass(wall.pscalar)
#    plot_effmass(point.pscalar)
#    plot_effmass(wall.vector)
#    plot_effmass(point.vector)
    
if __name__ == "__main__":
    main()
