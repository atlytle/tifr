import sys
import numpy as np
import itertools

from gammas import *
from read_HISQ import propagator_name, pion_correlator2, extract_t_fromfile
from read_HISQ import pion_correlator as pion_corr
from overlap_meson import spinmult

nx, ny, nz, nt = 24, 24, 24, 64
nc = 3
ns = 4

def HISQ_index(c1, c2, z, y, x):
    "Index of elements in (complexified, single timeslice) HISQ prop data."
    return c1*nz*ny*nx*nc + z*ny*nx*nc + y*nx*nc + x*nc + c2

def reshape_HISQ(propdata):
    "Reorder and add spin dofs."
    assert propdata.shape == (nx*ny*nz*nc*nc,)
    new_indices = [HISQ_index(c1,c2,z,y,x) 
                   for x,y,z,c1,c2 in itertools.product(range(nx), range(ny),
                                              range(nz), range(nc), range(nc))]
    # Reorder.
    tmp = propdata[new_indices]
    # Add spin dofs.
    tmp = np.kron(tmp, np.ones((ns*ns))).reshape(nx*ny*nz*nc*nc,ns,ns)
    return tmp

# Basic structure, details may be wrong!
def Omega(x,y,z,t):
    "Omega function used to wilsonize staggered propagator."
    tmp = id4
    if t % 2 == 0:
        tmp = np.dot(gt, tmp)
    if x % 2 == 0:
        tmp = np.dot(gx, tmp)
    if y % 2 == 0:
        tmp = np.dot(gy, tmp)
    if z % 2 == 0:
        tmp = np.dot(gz, tmp)
    return tmp
    
def Omega2(x,y,z,t):
    "Omega function used to wilsonize staggered propagator."
    tmp = id4
    if z % 2 == 1:
        tmp = np.dot(gz, tmp)
    if y % 2 == 1:
        tmp = np.dot(gy, tmp)
    if x % 2 == 1:
        tmp = np.dot(gx, tmp)
    if t % 2 == 1:
        tmp = np.dot(gt, tmp)
    return tmp
        
def Wilsonizer(x,y,z,t):
    "Factor to apply at each x,y,z,t,c1,c2 to convert staggered propagator."
    OmegaL = Omega2(x,y,z,t)
    # OmegaR = 1 so we do not include it here.
    return OmegaL

def wmatrix(t):
    return np.array([Wilsonizer(x,y,z,t) for x,y,z,c1,c2 in 
                    itertools.product(range(nx), range(ny), range(nz),
                    range(nc), range(nc))])

def pion_correlator(propname):
    '''HISQ pion correlator using Wilsonized propagator.

    This is mainly for illustrative purposes.  In particular the Wilson matrix
    cancels trivially in this case and can't be checked for correctness.
    '''
    tmp = extract_t_fromfile(propname, 0)  # Extract timeslice 0.
    tmp = reshape_HISQ(tmp)  # Shape is now (nx*ny*nz, ns, ns).
    tmp = wilson_matrix*tmp  # Multiply in spin.
    tmp1 = spinmult(g5, tmp)
    tmp1 = spinmult(g5, tmp1)
    assert (tmp == tmp1).all()  # For pion.
    tmp2 = np.conjugate(np.transpose(tmp, (0,2,1)))  # Transpose on spin.
    tmp2 = spinmult(g5, tmp2)
    tmp2 = spinmult(g5, tmp2)
    tmp2 = np.transpose(tmp2, (0,2,1))
    assert (np.conjugate(tmp2) == tmp1).all()  # For pion.
    return (tmp1*tmp2).astype(np.complex128).sum()

def kaon_correlator(prop1, prop2):
    '''Pseudoscalar correlator using two Wilsonized propagators.

    Mainly for illustrative purposes.
    '''
    
    # First propagator.
    tmp = extract_t_fromfile(prop1, 0)  # Extract timeslice 0.
    tmp = reshape_HISQ(tmp)  # Shape is now (nx*ny*nz, ns, ns).
    tmp = wilson_matrix*tmp  # Multiply in spin.
    tmp1 = spinmult(g5, tmp)
    tmp1 = spinmult(g5, tmp1)
    assert (tmp == tmp1).all()  # For pion.
    
    # Second propagator.
    tmp = extract_t_fromfile(prop2, 0)
    tmp = reshape_HISQ(tmp)
    tmp = wilson_matrix*tmp
    tmp2 = np.conjugate(np.transpose(tmp, (0,2,1)))  # Transpose on spin.
    tmp2 = spinmult(g5, tmp2)
    tmp2 = spinmult(g5, tmp2)
    tmp2 = np.transpose(tmp2, (0,2,1))
    
    return (tmp1*tmp2).astype(np.complex128).sum()

def meson_correlator(propname, g1, g2, t):
    '''General meson correlator using Wilsonized staggered propagator.

    Has the form Tr sum_x g2 g5 S^+ g5 g1 S, where S is the Wilsonized prop.
    '''
    # Construct Wilsonized propagator at timeslice t.
    tmp = extract_t_fromfile(propname, t)
    tmp = reshape_HISQ(tmp)
    tmp = wmatrix(t)*tmp

    tmp1 = spinmult(g1, tmp)
    tmp1 = spinmult(g5, tmp1)
    
    tmp2 = np.conjugate(np.transpose(tmp, (0,2,1)))  # Transpose on spin.
    tmp2 = spinmult(g5, tmp2)
    tmp2 = spinmult(g2, tmp2)
    tmp2 = np.transpose(tmp2, (0,2,1))
    
    print t, (tmp1*tmp2).astype(np.complex128).sum()

def main():
    # Load HISQ propagators.
    #prop1 = propagator_name('635', 1000)
    prop2 = propagator_name('0509', 1000)

    #print pion_correlator(prop2)
    print 'rho_y'
    for t in range(nt):
        meson_correlator(prop2, gy, gy, t)
    #print 'pion'
    #meson_correlator(prop2, g5, g5, 0)
    #meson_correlator(prop2, g5, g5, 1)


if __name__ == "__main__":
    sys.exit(main())
