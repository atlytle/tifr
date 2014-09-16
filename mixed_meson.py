import sys
import itertools
import numpy as np
from os.path import getsize

from gammas import *
from read_HISQ import (convert_to_complex, propagator_name,
                       pion_correlator2, extract_t_fromfile)
from hisq_meson import reshape_HISQ, wmatrix, wmatrixAM
from read_overlap import extract_t5
from read_overlap import prop_name as overlap_prop_loc
from overlap_meson import reshape_overlap, convert_to_aMILC, spinmult, smatrix
#from overlap_meson import spinmult

np.set_printoptions(precision=4)    
    
def mixed_pion_correlator(overlap_prop, HISQ_prop):
    "Pion correlator from one overlap propagator and one HISQ propagator."
    smat = smatrix()
    for t in range(nt):
        tmpO = extract_t5(overlap_prop, t)
        tmpO = reshape_overlap(tmpO)
        tmpO = convert_to_aMILC(tmpO)
        tmpO = smat*tmpO
        #tmpO = np.conjugate(np.transpose(tmpO, (0,2,1)))
        tmpH = extract_t_fromfile(HISQ_prop, t)
        tmpH = reshape_HISQ(tmpH)
        tmpH = wmatrixAM(t)*tmpH
        #tmpH = np.transpose(tmpH, (0,2,1))
        tmpH = np.conjugate(tmpH)
        print t, (tmpO*tmpH).astype(np.complex128).sum()

def mixed_meson_correlator(overlap_prop, HISQ_prop, g1, g2):
    '''General mixed meson correlator.
    '''
    smat = smatrix()
    for t in range(nt):
        tmpH = extract_t_fromfile(HISQ_prop, t)
        tmpH = reshape_HISQ(tmpH)
        tmpH = wmatrixAM(t)*tmpH
        tmpH = spinmult(g1,tmpH)
        tmpH = spinmult(g5M, tmpH)
        
        tmpO = extract_t5(overlap_prop, t)
        tmpO = reshape_overlap(tmpO)
        tmpO = convert_to_aMILC(tmpO)
        tmpO = smat*tmpO
        tmpO = np.conjugate(np.transpose(tmpO, (0,2,1)))
        tmpO = spinmult(g5M, tmpO)
        tmpO = spinmult(g2, tmpO)
        tmpO = np.transpose(tmpO, (0,2,1))
        print t, np.sum(tmpH*tmpO)

def mixed_meson_trace(overlap_prop, HISQ_prop, g1, g2):
    t = 0
    smat = smatrix()

    tmpH = extract_t_fromfile(HISQ_prop, t)
    tmpH = reshape_HISQ(tmpH)
    tmpH = wmatrixAM(t)*tmpH
    tmpH = spinmult(g1,tmpH)
    tmpH = spinmult(g5M, tmpH)
    tmp1 = tmpH[0]
    
    tmpO = extract_t5(overlap_prop, t)
    tmpO = reshape_overlap(tmpO)
    tmpO = convert_to_aMILC(tmpO)
    tmpO = smat*tmpO
    tmpO = np.conjugate(np.transpose(tmpO, (0,2,1)))
    tmpO = spinmult(g5M, tmpO)
    tmpO = spinmult(g2, tmpO)
    tmp2 = tmpO[0]
    tmpO = np.transpose(tmpO, (0,2,1))

    print tmpH[0].shape
    print tmpO[0].shape
    print np.trace(np.dot(tmp1,tmp2))
    print np.sum(tmpH[0]*tmpO[0])


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
    
    #mixed_pion_correlator(prop_ov2, prop2)  # Slow!
    #mixed_meson_correlator(prop_ov2, prop2, g5M, g5M)  # Slower!
    mixed_meson_trace(prop_ov2, prop2, g5M, g5M)

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
    

