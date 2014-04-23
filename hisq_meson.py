import sys
import numpy as np

from read_HISQ import propagator_name, pion_correlator2
from read_HISQ import pion_correlator as pion_corr
from read_mixed import (gx, gy, gz, gt, g5, id4,
                        reshape_HISQ, extract_t_fromfile, wilson_matrix)
from overlap_meson import spinmult

nx, ny, nz, nt = 24, 24, 24, 64
nc = 3
ns = 4

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

def main():
    # Load HISQ propagators.
    prop1 = propagator_name('635', 1000)
    prop2 = propagator_name('0509', 1000)

    print pion_correlator2(prop1, prop2)[0]  # Total agreement.
    print kaon_correlator(prop1, prop2)/4  # Total agreement.

if __name__ == "__main__":
    sys.exit(main())