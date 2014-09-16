import sys
import numpy as np
import itertools

from gammas import *
from read_overlap import extract_t5

nfloat = nx*ny*nz*nt*nc*nc*2  # Number of 4byte numbers expected.

def overlap_index(s1, c1, s2, c2, z, y, x):
    "Index of elements in (complexified, single timeslice) overlap prop data."
    return s1*nc*ns*nc*nz*ny*nx + c1*ns*nc*nz*ny*nx + s2*nc*nz*ny*nx +\
           c2*nz*ny*nx + z*ny*nx + y*nx + x

def reshape_overlap(propdata):
    "Reorder overlap data."
    assert propdata.shape == (ns*nc*ns*nc*nx*ny*nz,)
    new_indices = [overlap_index(s1,c1,s2,c2,z,y,x) for x,y,z,c1,c2,s1,s2 in
                      itertools.product(range(nx), range(ny), range(nz), 
                      range(nc), range(nc), range(ns), range(ns))]
    tmp = propdata[new_indices]
    return tmp.reshape((nx*ny*nz*nc*nc, ns, ns))

def spinmult(g, p):
    '''Multiply timeslice propagator by a gamma matrix.

    Preserves index structure.
    '''
    tmp = np.tensordot(g, p, axes=([1,1]))
    return np.transpose(tmp, (1,0,2))

def convert_to_aMILC(prop):
    '''Convert the Kentucky-convention overlap timeslice propagator
    to anti-MILC conventions.
    '''
    tmp = np.tensordot(T, prop, axes=([1,1]))
    tmp = np.transpose(tmp, (1,0,2))
    tmp = np.tensordot(tmp, hc(T), axes=([2,0]))
    return tmp

def smatrix():
    "To muliply off-diagonal elements by -1."
    a = np.array([[1,-1,-1,-1],
                  [-1,1,-1,-1],
                  [-1,-1,1,-1],
                  [-1,-1,-1,1]])
    smatrix = np.array([a for x,y,z,c1,c2 in
                         itertools.product(range(nx), range(ny), range(nz),
                         range(nc), range(nc))])
    return smatrix

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

    Has the form Tr sum_x g2 g5 S^+ g5 g1 S.
    '''
    corr = np.zeros((nt), dtype=complex)
    #smat = smatrix()
    for t in range(nt):
        tmp = extract_t5(propname, t)  # Extract timeslice t.
        tmp = reshape_overlap(tmp)  # Shape is now (nx*ny*nz, 4, 4).
        #tmp = smat*tmp
        tmp1 = spinmult(g1, tmp)
        tmp1 = spinmult(g5, tmp1)
        tmp2 = np.conjugate(np.transpose(tmp, (0,2,1)))  # Transpose on spin.
        tmp2 = spinmult(g5, tmp2)
        tmp2 = spinmult(g2, tmp2)
        tmp2 = np.transpose(tmp2, (0,2,1))
        corr[t] = np.sum(tmp1*tmp2)
        print t, corr[t]
    return corr



def main(filename):
    #print pion_correlator(filename)
    
    #print "pion"
    #pion = meson_correlator(filename, g5, g5)
    
    #print "rho_x"
    #rho_x = meson_correlator(filename, gx, gx)
    
    #print "rho_y"
    #rho_y = meson_correlator(filename, gy, gy)
    #print "rho_z"
    #rho_z = meson_correlator(filename, gz, gz)
    #rho = rho_x + rho_y + rho_z
    
    print "rho0"
    gzgt = np.dot(gz, gt)
    rho0 = meson_correlator(filename, gzgt, gzgt)
    print "pvy"
    g5gy = np.dot(g5, gy)
    pvy = meson_correlator(filename, g5gy, g5gy)
    print "B"
    gxgy = np.dot(gx, gy) 
    B = meson_correlator(filename, gxgy, gxgy)

if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))
