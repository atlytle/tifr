import numpy as np
from numpy import exp, cosh
from scipy.optimize import brentq, leastsq

from resample import JKsigma

# why not make these generic and just pass in fitfunc? also need initial guess.
def fit_double_cosh_osc(xarr, yarr, earr, T):
    "Fit (x,y,e) data to a double cosh form. Second cosh oscillates~(-1)^t."
    assert xarr.shape == yarr.shape == earr.shape
    
    fitfunc = lambda p, x: p[0]*cosh(p[1]*(x-T/2.)) + \
                          ((-1)**x)*p[2]*cosh(p[3]*(x-T/2.))
    errfunc = lambda p, x, y, err: (y-fitfunc(p,x))/err
    p0 = np.array([1.,1.,1.,1.])  # Initial guess.
    p1, success = leastsq(errfunc, p0, args=(xarr,yarr,earr), full_output=0)
    
    # Calculate chi^2.
    diffs = [errfunc(p1, x, y, err) for x, y, err in zip(xarr, yarr, earr)]
    chisq = np.sum(np.power(diffs,2))  # Not normalized.
    
    return p1, chisq
    
def fit_cfuns_double_cosh_osc(cfnc, ti, tf, T):
    xarr = np.array(range(T))
    earr = JKsigma(cfnc)  # These are held fixed throughout fitting.
    vals = [fit_double_cosh_osc(xarr[ti:tf], yarr[ti:tf], earr[ti:tf], T)
            for yarr in cfnc]  # Results of fit on each sample.
    chisq = vals[0][1]  # Chisq on central value fit.
    pvals = np.array([v[0] for v in vals])  # Parameter values of fits.
    
    return pvals[0], JKsigma(pvals), chisq


def fit_single_cosh_osc(xarr, yarr, earr, T):
    "Fit (x, y, e) data to (A+(-B)^t) cosh form."
    assert xarr.shape == yarr.shape == earr.shape
    
    fitfunc = lambda p, x: (p[0]+((-1)**x)*p[1])*cosh(p[2]*(x-T/2.))
    errfunc = lambda p, x, y, err: (y-fitfunc(p,x))/err
    p0 = np.array([1.,1.,1.])  # Initial guess.
    p1, success = leastsq(errfunc, p0, args=(xarr,yarr,earr), full_output=0)
    
    # Calculate chi^2.
    diffs = [errfunc(p1, x, y, err) for x, y, err in zip(xarr, yarr, earr)]
    chisq = np.sum(np.power(diffs,2))  # Not normalized.
    
    return p1, chisq

def fit_cfuns(cfnc, ti, tf, T, fitfunc):
    xarr = np.array(range(T))
    earr = JKsigma(cfnc)
    vals = [fitfunc(xarr[ti:tf], yarr[ti:tf], earr[ti:tf], T)
            for yarr in cfnc]  # Results of fit on each sample.
    chisq = vals[0][1]  # Chisq on central value fit.
    pvals = np.array([v[0] for v in vals])  # Parameter values of fits.
    
    return pvals[0], JKsigma(pvals), chisq

def fit_cfuns_through(cfnc, ti, tf, T, fitfunc):
    """Version of fit_cfuns that returns all the resampled results."""
    xarr = np.array(range(T))
    earr = JKsigma(cfnc)
    vals = [fitfunc(xarr[ti:tf], yarr[ti:tf], earr[ti:tf], T)
            for yarr in cfnc]  # Results of fit on each sample.
    chisq = vals[0][1]  # Chisq on central value fit.
    pvals = np.array([v[0] for v in vals])  # Parameter values of fits.
    
    return pvals, chisq

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
    
def naive_effmass(cfnc, foldQ=False):
    '''Compute naive effective mass of correlation functions.'''
    if foldQ:
        cfnc = fold(cfnc)
    cplus = np.roll(cfnc, 1, axis=1).real
    return np.abs(np.log(cfnc/cplus))
