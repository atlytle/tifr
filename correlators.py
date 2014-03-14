import numpy as np
from numpy import exp, cosh
from scipy.optimize import brentq, leastsq

from resample import JKsigma

# why not make these generic and just pass in fitfunc?
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
