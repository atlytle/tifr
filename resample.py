import numpy as np

from random import randrange

'''
Resampling distributions.

It would be nice to make these datatype agnostic, but numpy arrays
don't have append() while lists don't have vstack().  Thus the current
examples are targeted to (2d) numpy arrays.

'''

def yield_JK_sample_np(arr):
    Nconf = arr.shape[0]
    orig = arr.copy()
    for i in range(Nconf):
        x, sample = arr[0:1], arr[1:]
        yield sample
        arr = np.vstack((sample, x))
    assert (arr == orig).all()

def yield_bootstrap_sample_np(arr, Nboot):
    Nconf = arr.shape[0]
    for n in range(Nboot):
        indices = [randrange(Nconf) for i in range(Nconf)]
        sample = [arr[i] for i in indices]
        yield np.array(sample)
        
