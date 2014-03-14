import numpy as np

from random import randrange

'''
Resampling distributions.

It would be nice to make these datatype agnostic, but numpy arrays
don't have append() while lists don't have vstack().  Thus the current
examples are targeted to (2d) numpy arrays.

'''

def yield_JK_sample_np(arr):
    "Yield jackknife samples from rows of a numpy array."
    Nconf = arr.shape[0]
    orig = arr.copy()
    for i in range(Nconf):
        x, sample = arr[0:1], arr[1:]
        yield sample
        arr = np.vstack((sample, x))
    assert (arr == orig).all()
    
def JK_block(arr):
    "Average and jackknife averages of rows of a numpy array."
    ave = lambda x: np.average(x, axis=0)
    return np.vstack((ave(arr), map(ave, yield_JK_sample_np(arr))))
    
def JKsigma(cfnc):
    '''Errorbars given correlator average and jackknife correlators.'''
    ave, JKvals = cfnc[0], cfnc[1:]
    N = len(JKvals)
    diffs = (ave-JKvals)*(ave-JKvals)
    return np.sqrt(np.sum(diffs, axis=0)*(1-1./N))

def yield_bootstrap_sample_np(arr, Nboot):
    "Yield Nboot bootstrap samples from rows of a numpy array."
    Nconf = arr.shape[0]
    for n in range(Nboot):
        indices = [randrange(Nconf) for i in range(Nconf)]
        sample = [arr[i] for i in indices]
        yield np.array(sample)
        
def block(meas, sample, arr):
    '''Generalized version for resampled measurements.
    
    e.g.  block(lambda x: np.average(x, axis=0), yield_JK_sample_np, arr)
          == JK_block(arr).  Still assumes numpy datatypes though.
          
    '''
    return np.vstack((meas(arr), map(meas, sample(arr))))
        
