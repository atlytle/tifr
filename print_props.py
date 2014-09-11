"""Manual printing of propagators for debugging/checking purposes."""

import sys
import itertools
import numpy as np
np.set_printoptions(precision=4)

#from os.path import getsize

from read_HISQ import propagator_name

#from read_overlap import extract_t5
from read_mixed import (extract_t_fromfile, HISQ_index, reshape_HISQ, 
                        wmatrix, T)
                        #gx, gy, gz, gt, g5, id4)

nx, ny, nz, nt = 24, 24, 24, 64
nc = 3
ns = 4
nfloat = nx*ny*nz*nt*nc*nc*2  # Number of 4byte numbers expected.

def main(argv=None):                                                                                          
    hprop = propagator_name('0509', 1000)
    for x,y,z,t in (0,0,0,0), (1,0,0,0), (0,1,0,0), (0,0,1,0), (0,0,0,1),\
                   (1,1,0,0), (1,0,1,0), (1,0,0,1), (0,1,1,0),\
                   (0,1,0,1), (0,0,1,1), (1,1,1,0), (1,1,0,1),\
                   (1,0,1,1), (0,1,1,1), (1,1,1,1):
        print 'x={0}, y={1}, z={2}, t={3}'.format(x,y,z,t)
        print '--------------------------'
        tmpH = extract_t_fromfile(hprop, t)
        tmpH2 = reshape_HISQ(tmpH)
        tmpH3 = wmatrix(t)*tmpH2
        for c1 in range(nc):
            for c2 in range(nc):
                print 'c1={0}, c2={1}'.format(c1,c2)
                print '--'
                id = x*ny*nz*nc*nc + y*nz*nc*nc + z*nc*nc + c1*nc + c2
                print 'Original data:  {0:.4f}'.format(tmpH2[id,0,0])
                print 'Wilsonized data (MILC):'
                for s1 in range(ns):
                    for s2 in range(ns):
                        print 's1={0}, s2={1}:  {2:.4f}'.format(s1, s2, 
                                                            tmpH3[id,s1,s2])
                print ''
        print '\n'

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
