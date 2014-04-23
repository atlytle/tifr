import sys
import numpy as np
from read_overlap import extract_t5
from read_mixed import gx, gy, gz, gt, g5, id4, reshape_overlap

nx, ny, nz, nt = 24, 24, 24, 64
nc = 3
ns = 4
nfloat = nx*ny*nz*nt*nc*nc*2  # Number of 4byte numbers expected.

def spinmult(g, p):
    '''Multiply timeslice propagator by a gamma matrix.

    Preserves index structure.
    '''
    tmp = np.tensordot(g, p, axes=([1,1]))
    return np.transpose(tmp, (1,0,2))

def pion_correlator(propname):
    '''Overlap pion correlator using spin dofs.

    This is mainly for illustrative purposes and to check other cases.
    '''
    tmp = extract_t5(propname, 0)  # Extract timeslice 0.
    tmp = reshape_overlap(tmp)  # Shape is now (nx*ny*nz, 4, 4).
    tmp1 = spinmult(g5, tmp)
    tmp1 = spinmult(g5, tmp1)
    assert (tmp == tmp1).all()  # For pion.
    tmp2 = np.conjugate(np.transpose(tmp, (0,2,1)))  # Transpose on spin.
    tmp2 = spinmult(g5, tmp2)
    tmp2 = spinmult(g5, tmp2)
    tmp2 = np.transpose(tmp2, (0,2,1))
    assert (np.conjugate(tmp2) == tmp1).all()  # For pion.
    return np.sum(tmp1*tmp2)

def meson_correlator(propname, g1, g2):
    '''General meson correlator.

    Has the form Tr sum_x g1 g5 S^+ g5 g2 S.
    '''
    for t in range(nt):
        tmp = extract_t5(propname, t)  # Extract timeslice t.
        tmp = reshape_overlap(tmp)  # Shape is now (nx*ny*nz, 4, 4).
        tmp1 = spinmult(g5, tmp)
        tmp1 = spinmult(g1, tmp1)
        tmp2 = np.conjugate(np.transpose(tmp, (0,2,1)))  # Transpose on spin.
        tmp2 = spinmult(g5, tmp2)
        tmp2 = spinmult(g2, tmp2)
        tmp2 = np.transpose(tmp2, (0,2,1))
        print t, np.sum(tmp1*tmp2)



def main(filename):
    print pion_correlator(filename)
    meson_correlator(filename, gx, gx)

if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))