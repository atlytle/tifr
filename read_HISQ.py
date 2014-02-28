import sys
import struct
import numpy as np

nx, ny, nz, nt = 24, 24, 24, 64
nc = 3
                
def extract_t(clist, t):
# c c t z y x
    result = []
    for c in range(3):
        i = c*nt*nz*ny*nx*nc + t*nz*ny*nx*nc
        result += clist[i:i+nz*ny*nx*nc]
    assert len(result) == nc*nc*nz*ny*nx
    return result
        
                
def main(filename):
    num_float = nx*ny*nz*nt*nc*nc*2  # Number of 4byte numbers expected.
    data = np.fromfile(filename, dtype='<f')
    data = data[28:num_float+28]  # Discard header.
    assert len(data) == num_float
    print data[:10]
    re = data[::2]
    im = data[1::2]
    c = map(complex, re, im)
    assert len(c) == nx*ny*nz*nt*nc*nc
#    for t in range(nt):
#        tmp = extract_t(c, t)
#        tmp = np.array(tmp)
#        if t==0:
#            i = 0
#            for x in tmp:
#                i += 1
#                if abs(x) > 1: 
#                    print x
#                    print i
#        print t, sum(np.square(np.abs(tmp)))
         
if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))

