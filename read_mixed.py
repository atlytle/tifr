import sys
import itertools
import numpy as np

from os.path import getsize

from read_HISQ import (convert_to_complex, propagator_name,
                       pion_correlator2, extract_t)
from read_overlap import extract_t5
#from overlap_meson import spinmult

nx, ny, nz, nt = 24, 24, 24, 64
nc = 3
ns = 4
nfloat = nx*ny*nz*nt*nc*nc*2  # Number of 4byte numbers expected.

# Define gamma matrices.  Need to CHECK overlap & hisq conventions.

def hc(M):
    "Hermitian conjugate."
    return np.conjugate(np.transpose(M))  # Transpose on spin.


# Kentucky format gamma matrices.
# The overlap propagators are written in this convention.

gx = np.array([[0, 0, 0, -1j],
               [0, 0, -1j, 0],
               [0, 1j, 0, 0],
               [1j, 0, 0, 0]])

gy = -np.array([[0, 0, 0, 1],
               [0, 0, -1, 0],
               [0, -1, 0, 0],
               [1, 0, 0, 0]])

gz = np.array([[0, 0, -1j, 0],
               [0, 0, 0, 1j],
               [1j, 0, 0, 0],
               [0, -1j, 0, 0]])

gt = np.array([[1, 0, 0, 0],
               [0, 1, 0, 0],
               [0, 0, -1, 0],
               [0, 0, 0, -1]])
               
g5 = np.array([[0, 0, 1, 0],
               [0, 0, 0, 1],
               [1, 0, 0, 0],
               [0, 1, 0, 0]])

id4 = np.identity(4, dtype=complex)
id3 = np.identity(3, dtype=complex)

# Transformation from Kentucky to "Wilson"
T = np.array([[0, -1, 0, 1],
              [1, 0, -1, 0],
              [0, 1, 0, 1],
              [-1, 0, -1, 0]])/np.sqrt(2)

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
        
omega_matrix = np.array([Omega(x,y,z,0) for x,y,z,c1,c2 in 
                          itertools.product(range(nx), range(ny), range(nz),
                          range(nc), range(nc))])
        
def Wilsonizer(x,y,z,t):
    "Factor to apply at each x,y,z,t,c1,c2 to convert staggered propagator."
    OmegaL = reduce(np.dot, [id4, Omega2(x,y,z,t), id4])  # hc[T] == -T.
    #OmegaR = hc(Omega2(0,0,0,0))  # This is very inefficient...unnecessary.
    #return reduce(np.dot, [OmegaL, id4, OmegaR])
    return OmegaL

def wmatrix(t):
    return np.array([Wilsonizer(x,y,z,t) for x,y,z,c1,c2 in 
                    itertools.product(range(nx), range(ny), range(nz),
                    range(nc), range(nc))])

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
    
def overlap_prop_loc(m, config):
    root = '/user2/atlytle/overlap/L24T64/'
    return root + 'l2464f211b600m0102m0509m635a.'\
    '{0}.gf_hyp_prop_src_wall_t00_m{1:7f}'.format(config, m)
    
def HISQ_index(c1, c2, z, y, x):
    "Index of elements in (complexified, single timeslice) HISQ prop data."
    return c1*nz*ny*nx*nc + z*ny*nx*nc + y*nx*nc + x*nc + c2
    
def overlap_index(s1, c1, s2, c2, z, y, x):
    "Index of elements in (complexified, single timeslice) overlap prop data."
    return s1*nc*ns*nc*nz*ny*nx + c1*ns*nc*nz*ny*nx + s2*nc*nz*ny*nx +\
           c2*nz*ny*nx + z*ny*nx + y*nx + x
    
def extract_t_fromfile(filename, t):
    '''Extract timeslice t from HISQ propagator.
    
    Return complex values ordered like c1 z y x c2.'''
    
    # Check file size.
    if getsize(filename) != (nfloat+31)*4:
        raise Exception('{0} does not have the expected size.'.format(filename))
    
    f = open(filename, 'rb')
    
    # Read header.  This consists of 22 4byte values.
    data = np.fromfile(f, dtype='<i', count=5)  # 5 floats.
    timestamp = np.fromfile(f, dtype='<c', count=64)  # 16 floats.
    order = np.fromfile(f, dtype='<i', count=1)  # 1 float.
    

    
    # Extract data at timeslice t for each color chunk.
    result = np.array([])
    for c in range(nc):
        color = np.fromfile(f, dtype='<i', count=1)
        some_nums = np.fromfile(f, dtype='<f', count=2) 
        data = np.fromfile(f, dtype='<f', count=nfloat/3)
        tmp = extract_t(data, t) 
        result = np.concatenate((result,tmp))
        
    f.close()
    return convert_to_complex(result)
    
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
    
def reshape_overlap(propdata):
    "Reorder overlap data."
    assert propdata.shape == (ns*nc*ns*nc*nx*ny*nz,)
    new_indices = [overlap_index(s1,c1,s2,c2,z,y,x) for x,y,z,c1,c2,s1,s2 in
                      itertools.product(range(nx), range(ny), range(nz), 
                      range(nc), range(nc), range(ns), range(ns))]
    tmp = propdata[new_indices]
    return tmp.reshape((nx*ny*nz*nc*nc, ns, ns))
    
    
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
    

