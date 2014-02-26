import sys
import struct
import numpy as np

nx, ny, nz, nt = 24, 24, 24, 64
nc = 3

def yield_bytes(filename):
    with open(filename, "rb") as f:
        while True:
            byte = f.read(4)
            if byte:
                yield byte
            else:
                break
                
def yield_bytes2(filename):
    with open(filename, "rb") as f:
        while True:
            byte = f.read(4)
            if byte:
                yield struct.unpack_from('f', byte)[0]
            else:
                break
                
def extract_t(clist, t):
# c c t z y x
    result = []
    for c in range(3):
        i = c*nt*nz*ny*nx*nc + t*nz*ny*nx*nc
        result += clist[i:i+nz*ny*nx*nc]
    assert len(result) == nc*nc*nz*ny*nx
    return result
        
                
def main(filename):
    nx, ny, nz, nt = 24, 24, 24, 64
    nc = 3
    num_float = nx*ny*nz*nt*nc*nc*2  # Number of 4byte numbers expected.
    floats = [b for b in yield_bytes2(filename)]
    floats = floats[28:num_float+28]  # Discard any header.
    print floats[:50]
    #print floats[-50:]
    assert len(floats) == num_float
    re = floats[::2]
    im = floats[1::2]
    c = map(complex, re, im)
    assert len(c) == nx*ny*nz*nt*nc*nc
    #print c[0:3]
    for t in range(nt):
        tmp = extract_t(c, t)
        tmp = np.array(tmp)
#        if t==0:
#            i = 0
#            for x in tmp:
#                i += 1
#                if abs(x) > 1: 
#                    print x
#                    print i
        print t, sum(np.square(np.abs(tmp)))
         
if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))

