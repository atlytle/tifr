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
    
def plot_effmass(cfnc):
    pass
    
def naive_effmass(cfnc):
    '''Compute effective mass of correlation function.'''
    pass
    
    
def main():
    ms = 0.0495
    wall = OverlapWall(0.55, ms)
    point = OverlapPoint(0.551, ms)
    
    plot_correlator(wall.pscalar)
    plot_correlator(point.pscalar)
    plot_correlator(wall.vector)
    plot_correlator(point.vector)

if __name__ == "__main__":
    main()
