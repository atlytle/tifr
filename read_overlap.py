import sys
import struct
import itertools
import numpy as np

nx, ny, nz, nt = 32, 32, 32, 96
V = nx*ny*nz  # Spatial volume.
nc = 3  # N_c.
ns = 4  # N_spin.
ndouble = V*nt*nc*nc*ns*ns*2  # Number of 8byte numbers expected.

def extract_t(data, t):
    block = V*nc*nc*ns*ns*2
    return data[t*block:t*block + block]  

def reorder(data):
    "Reorder the data in a sensible way."
    num = -1
    assert len(data) == ndouble
    result = np.zeros((ndouble))
    loop = itertools.product(range(nt), range(ns), range(nc), range(2), 
                             range(ns), range(nc), range(nz), range(ny), 
                             range(nx))
    for t, s2, c2, ri, s1, c1, z, y, x in loop:
        num += 1
        old = x + nx*y + nx*ny*z + nx*ny*nz*c1 + nx*ny*nz*nc*s1 + \
              nx*ny*nz*nc*ns*ri + nx*ny*nz*nc*ns*2*c2 + \
              nx*ny*nz*nc*ns*2*nc*s2 + nx*ny*nz*nc*ns*2*nc*ns*t
        assert old == num
        new = ri + 2*c2 + 2*nc*s2 + 2*nc*ns*c1 + 2*nc*ns*nc*s1 + \
              2*nc*ns*nc*ns*x + 2*nc*ns*nc*ns*nx*y + 2*nc*ns*nc*ns*nx*ny*z +\
              2*nc*ns*nc*ns*nx*ny*nz*t
        result[new] = data[old]
    assert num == ndouble - 1
    assert 0 not in result
    return result

                
def main(filename):
    for t in range(nt):
        result = []
        for i in range(2*nc*nc*ns*ns):
            with open(filename, "rb") as f:
                f.seek(i*8*nt*V+8*t*V, 0)
                data = np.fromfile(f, dtype='>d', count=V)
                result.append(np.sum(np.square(data)))
            
        print t, sum(result)
    
         
if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))

