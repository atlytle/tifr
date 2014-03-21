import sys
import itertools
import numpy as np

from os.path import getsize

from read_HISQ import convert_to_complex

nx, ny, nz, nt = 24, 24, 24, 64
nc = 3
nfloat = nx*ny*nz*nt*nc*nc*2  # Number of 4byte numbers expected.

def propagator_name(m, config):
    root = '/user2/atlytle/staggered'
    return root + '/hisq_{0}/prop_bin_hisq_{0}.{1}'.format(m, config)
    
def reorder_HISQ(c1, c2, z, y, x):
    return c1*nz*ny*nx*nc + z*ny*nx*nc + y*nx*nc + x*nc + c2
    
def extract_t(data, t):
    "Extract chunk of data corresponding to timeslice t."
    # Loop structure: t z y x c r/i
    result = data[t*nz*ny*nx*nc*2:t*nz*ny*nx*nc*2 + nz*ny*nx*nc*2]
    assert len(result) == nz*ny*nx*nc*2
    return result
    
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
        print tmp.shape   
        result = np.concatenate((result,tmp))
        print result.shape
        
    f.close()
    return convert_to_complex(result)


def main(argv=None):                                              
    new_indices = [reorder_HISQ(c1,c2,z,y,x) 
                   for c1, c2, z, y, x in itertools.product(range(nc), range(nc),
                                              range(nz), range(ny), range(nx))]
                                              
    print new_indices[0:10]
    
    propname = propagator_name('0509', 1000)
    extract_t_fromfile(propname, 0)

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
    

