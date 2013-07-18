import optparse
import sys
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

def plot_correlator(cfnc, save=False, name='', title=None):
    ave, sigma = cfnc[0], JKsigma(cfnc)
    p.figure()
    if title is not None:
        p.title(title)
    p.xlabel('$t$')
    p.ylabel('$C[t]$')
    p.yscale('log')
    p.errorbar(range(len(ave)), ave.real,sigma.real, fmt='k-')
    if save:
        p.savefig(name)
    else:
        p.show()
    
def fit_twopoint(xarr, yarr, earr, T=96.): #might add t1, t2 here
    '''Fit two-point correlator to the cosh form.'''
    assert xarr.shape == yarr.shape == earr.shape
    
    fitfunc = lambda p, x: p[0]*(2*exp(-p[1]*T/2.)*cosh(p[1]*(x-T/2.)))
    errfunc = lambda p, x, y, err: (y - fitfunc(p,x))/err
    
    p0 = np.array([0.5, 0.5])  # Initial guess.
    p1, success = leastsq(errfunc, p0, args=(xarr,yarr,earr),
                                   full_output=0)
    # Calculate chi^2.
    diffs = [errfunc(p1, x, y, err) for x, y, err in zip(xarr, yarr, earr)]
    chisq = np.sum(np.power(diffs,2))  # Not normalized.
    
    return p1, chisq

def fit_twopoint_cfuns(cfnc, ti, tf, T):
    '''Fit two-point correlators to cosh form from t_initial to t_final. 
    
    Returns (A, m), (sig_A, sig_m), chisq.
    '''
    xarr = np.array(range(T))
    earr = JKsigma(cfnc)  # These are held fixed throughout fitting.
    vals = [fit_twopoint(xarr[ti:tf], yarr[ti:tf], earr[ti:tf], T)
            for yarr in cfnc]
    chisq = vals[0][1]  # Chisq on central value fit.
    vals = [v[0] for v in vals]  # Parameter values of fits.
    
    return vals[0], JKsigma(np.array(vals)), chisq
    
def fit_twopoint_cfunsJK(cfnc, ti, tf, T):
    '''Fit two-point correlators to cosh form from t_initial to t_final. 
    
    More resampling friendly version:
    Returns [(A, m)] for total-sample (0th element) and JK samples.
    '''
    xarr = np.array(range(T))
    earr = JKsigma(cfnc)  # These are held fixed throughout fitting.
    vals = [fit_twopoint(xarr[ti:tf], yarr[ti:tf], earr[ti:tf], T)
            for yarr in cfnc]
    vals = [v[0] for v in vals]  # Parameter values of fits.
    
    return vals
    
    
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
    
def plot_effmass2(cfnc, save=False, name=''):
    effmass = cosh_effmass(cfnc)
    ave, sigma = effmass[0], JKsigma(effmass)
    p.figure()
    p.xlabel('$t$')
    p.ylabel('$|\\frac{C[t+1]}{C[t]}|$')
    print len(ave)
    p.errorbar(range(len(ave)), ave, sigma, fmt='ko') # Fix. t-axis misaligned?
    if save:
        p.savefig(name)
    else:
        p.show()
    
def plot_effmass(cfnc, save=False, name='', title=None, yrange=[]):
    effmass = naive_effmass(cfnc)
    ave, sigma = effmass[0], JKsigma(effmass)
    p.figure()
    if title is not None:
        p.title(title)
    if yrange:
        p.ylim(yrange)
    p.xlabel('$t$')
    p.ylabel('$|\\frac{C[t+1]}{C[t]}|$')
    p.errorbar(range(len(ave)), ave, sigma, fmt='ko')
    if save:
        p.savefig(name)
    else:
        p.show()
        
def parse_args(argv):
    parser = optparse.OptionParser()
    parser.add_option('-p', '--plot', action='store_true', dest='plot',
                      help = 'Plot results.')
    parser.add_option('-s', '--save', action='store_true', dest='save',
                      help = 'Save the plots.')
    options, args = parser.parse_args(argv)
    return options
    
def main(argv=None):
    options = parse_args(argv)
    root = '/Users/atlytle/Dropbox/TIFR/figs/'
    
    ms = 0.0495
    wall = OverlapWall(0.55, ms)
    point = OverlapPoint(0.551, ms)

    p1, err, chisq = fit_twopoint_cfuns(wall.pscalar.real, 10, 40, wall.T)
    print 'Wall-Point pseudoscalar:'
    print 'A=', p1[0], '+/-', err[0]
    print 'm=', p1[1], '+/-', err[1]
    print 'chisq=', chisq, '\n'
    
    p1, err, chisq = fit_twopoint_cfuns(point.pscalar.real, 10, 40, wall.T)
    print 'Point-Point pseudoscalar:'    
    print 'A=', p1[0], '+/-', err[0]
    print 'm=', p1[1], '+/-', err[1]
    print 'chisq=', chisq, '\n'
    
    p1, err, chisq = fit_twopoint_cfuns(wall.vector.real, 10, 40, wall.T)
    print 'Wall-Point vector:'
    print 'A=', p1[0], '+/-', err[0]
    print 'm=', p1[1], '+/-', err[1]
    print 'chisq=', chisq, '\n'
    
    p1, err, chisq = fit_twopoint_cfuns(point.vector.real, 10, 40, wall.T)
    print 'Point-Point vector:'    
    print 'A=', p1[0], '+/-', err[0]
    print 'm=', p1[1], '+/-', err[1]
    print 'chisq=', chisq, '\n'
    
    p1, err, chisq = fit_twopoint_cfuns(-wall.a4a4.real, 10, 40, wall.T)
    print 'Wall-Point a4-a4:'
    print 'A=', p1[0], '+/-', err[0]
    print 'm=', p1[1], '+/-', err[1]
    print 'chisq=', chisq, '\n'
    
    p1, err, chisq = fit_twopoint_cfuns(-point.a4a4.real, 10, 40, wall.T)
    print 'Point-Point a4-a4:'    
    print 'A=', p1[0], '+/-', err[0]
    print 'm=', p1[1], '+/-', err[1]
    print 'chisq=', chisq, '\n'
    
    
#    plot_correlator(wall.pscalar, options.save, root+'wall_pscalar_corr.pdf')
#    plot_correlator(point.pscalar, options.save, root+'point_pscalar_corr.pdf')
#    plot_correlator(wall.vector, options.save, root+'wall_vector_corr.pdf')
#    plot_correlator(point.vector, options.save, root+'point_vector_corr.pdf')
#    plot_correlator(-point.a4a4, options.save, root+'point_a4a4_corr.pdf')
#    plot_correlator(-wall.a4a4, options.save, root+'wall_a4a4_corr.pdf')
#    plot_correlator(point.psa4, options.save, root+'point_psa4_corr.pdf')
#    plot_correlator(wall.psa4, options.save, root+'wall_psa4_corr.pdf')
#    plot_correlator(point.a4ps, options.save, root+'point_a4ps_corr.pdf')
#    plot_correlator(wall.a4ps, options.save, root+'wall_a4ps_corr.pdf')
#    plot_effmass(wall.pscalar, options.save, root+'wall_pscalar_meff_naive.pdf')
#    #plot_effmass2(wall.pscalar)
#    plot_effmass(point.pscalar, options.save, root+'point_pscalar_meff_naive.pdf')
#    plot_effmass(wall.vector, options.save, root+'wall_vector_meff_naive.pdf')
#    plot_effmass(point.vector, options.save, root+'point_vector_meff_naive.pdf')
#    plot_effmass(wall.a4a4, options.save, root+'wall_a4a4_meff_naive.pdf')
#    plot_effmass(point.a4a4, options.save, root+'point_a4a4_meff_naive.pdf')
#    plot_effmass(wall.psa4, options.save, root+'wall_psa4_meff_naive.pdf')
#    plot_effmass(point.psa4, options.save, root+'point_psa4_meff_naive.pdf')
#    plot_effmass(wall.a4ps, options.save, root+'wall_a4ps_meff_naive.pdf')
#    plot_effmass(point.a4ps, options.save, root+'point_a4ps_meff_naive.pdf')  
  
    return 0
    
if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
