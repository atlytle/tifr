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
    
def extract_t2(clist, t):
# c t z y x
    result = []
    #print clist.shape
    for c in range(3):
        i = c*nt*nz*ny*nx*2 + t*nz*ny*nx*2
        result =  np.hstack((result,clist[i:i+nz*ny*nx]))
    assert len(result) == nc*nz*ny*nx
    return result
    
def extract_t3(clist, t):
# c t z y x
    result = clist[t*nz*ny*nx*nc*2:t*nz*ny*nx*nc*2 + nz*ny*nx*nc*2]
    assert len(result) == nc*nz*ny*nx*2
    return result
        
                
def main(filename):
    nfloat = nx*ny*nz*nt*nc*nc*2  # Number of 4byte numbers expected.

    with open(filename, "rb") as f:
        # Read header.  This consists of 22 4byte values.
        data = np.fromfile(f, dtype='<i', count=5)  # 5 floats.
        timestamp = np.fromfile(f, dtype='<c', count=64)  # 16 floats.
        order = np.fromfile(f, dtype='<i', count=1)  # 1 float.
        print data
        print timestamp
        print order
        # Read data and calculate pion correlator.
        # Each chunk of data is separated by an int specifying source color,
        # and two other 4byte values.  Thus there are 31 'extra' 4bytes.
        result = np.zeros((nt))
        for c in range(nc):
            print 'color', np.fromfile(f, dtype='<i', count=1)
            print 'some nums', np.fromfile(f, dtype='<f', count=2) 
            x = np.fromfile(f, dtype='<f', count=nfloat/3)
            # Add intermediate results to correlator.
            # The sum is done in double precision.
            for t in range(nt):
                result[t] += np.square(extract_t3(x, t)).astype(np.float64).sum()
        for t in range(nt):
            print t, result[t]
         
if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))

