import sys
import itertools
import numpy as np
from os.path import getsize

from gammas import *
from read_HISQ import (convert_to_complex, propagator_name,
                       pion_correlator2, extract_t_fromfile)
from hisq_meson import reshape_HISQ, wmatrix
from read_overlap import extract_t5
from read_overlap import prop_name as overlap_prop_loc
from overlap_meson import reshape_overlap
#from overlap_meson import spinmult

nx, ny, nz, nt = 24, 24, 24, 64
nc = 3
ns = 4
nfloat = nx*ny*nz*nt*nc*nc*2  # Number of 4byte numbers expected.

def hc(M):
    "Hermitian conjugate."
    return np.conjugate(np.transpose(M)) 

def odmatrix():
    "Muliply off-diagonal elements by -1."
    a = np.array([[1,-1,-1,-1],
                  [-1,1,-1,-1],
                  [-1,-1,1,-1],
                  [-1,-1,-1,1]])
    odmatrix = np.array([a for x,y,z,c1,c2 in
                         itertools.product(range(nx), range(ny), range(nz),
                         range(nc), range(nc))])
    return odmatrix
    
    
def mixed_pion_correlator(overlap_prop, HISQ_prop):
    "Pion correlator from one overlap propagator and one HISQ propagator."
    for t in range(nt):
        tmpO = extract_t5(overlap_prop, t)
        tmpO = reshape_overlap(tmpO)
        tmpO = odmatrix()*tmpO
        tmpH = extract_t_fromfile(HISQ_prop, t)
        tmpH = reshape_HISQ(tmpH)
        tmpH = wmatrix(t)*tmpH
        print t, (tmpO*np.conj(tmpH)).astype(np.complex128).sum()

def mixed_meson_correlator(overlap_prop, HISQ_prop):
    "Unfinished."
    for t in range(nt):
        tmpH = extract_t_fromfile(HISQ_prop, t)
        tmpH = reshape_HISQ(tmp)
        tmpO = extract_t5(overlap_prop, t)
        tmpO = reshape_overlap(tmpO) 

def main(argv=None):                                                                                          
    
    # Load overlap propagators.
    #prop_ov1 = overlap_prop_loc(0.038, 1000)
    prop_ov2 = overlap_prop_loc(0.0731, 1000)
    
    # Load HISQ propagators.
    #prop1 = propagator_name('635', 1000)
    prop2 = propagator_name('0509', 1000)

    #t = 5
    #tmpO = extract_t5(prop_ov2, t)
    #tmpO = reshape_overlap(tmpO)
    #tmpH = extract_t_fromfile(HISQ_prop, t)
    #tmpH = reshape_HISQ(tmpH)
    #tmpH = wilson_matrix*tmpH
    
    mixed_pion_correlator(prop_ov2, prop2)  # Slow!

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
    

