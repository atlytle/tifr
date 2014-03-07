import sys
import struct
import numpy as np

from os.path import getsize

nx, ny, nz, nt = 24, 24, 24, 64
nc = 3
nfloat = nx*ny*nz*nt*nc*nc*2  # Number of 4byte numbers expected.

def correlator_name(m):
    return 'HHpion_l2464_m{0}_m{0}.npy'.format(m)
    
def extract_t(data, t):
    "Extract chunk of data corresponding to timeslice t."
    # Loop structure: t z y x c r/i
    result = data[t*nz*ny*nx*nc*2:t*nz*ny*nx*nc*2 + nz*ny*nx*nc*2]
    assert len(result) == nc*nz*ny*nx*2
    return result

def pion_correlator(filename):
    "Construct pion correlator from propagator."
    
    # Check file size.
    if getsize(filename) != (nfloat+31)*4:
        raise Exception('{0} does not have the expected size.'.format(filename))
    correlator = np.zeros((nt))
    
    f = open(filename, 'rb')
    
    # Read header.  This consists of 22 4byte values.
    data = np.fromfile(f, dtype='<i', count=5)  # 5 floats.
    timestamp = np.fromfile(f, dtype='<c', count=64)  # 16 floats.
    order = np.fromfile(f, dtype='<i', count=1)  # 1 float.
    
    # Read data and calculate pion correlator.
    # Each chunk of data is separated by an int specifying source color,
    # and two other 4byte values.  Thus there are 31 'extra' 4bytes.
    result = np.zeros((nt))
    for c in range(nc):
        color = np.fromfile(f, dtype='<i', count=1)
        some_nums = np.fromfile(f, dtype='<f', count=2) 
        data = np.fromfile(f, dtype='<f', count=nfloat/3)
        # Add intermediate results to correlator.
        # The sum is done in double precision.
        for t in range(nt):
            correlator[t] += np.square(extract_t(data, t)).astype(np.float64).sum()
    
    f.close()
    
    return correlator
            
                
def main(files):
    # Basic check on inputs.
    head0, config0 = files[0].split('.')
    mspec = head0.split('_')[-1]  # Mass specifier.
    for f in files:
        head, config = f.split('.')
        if head != head0:
            print "You might not want to combine these!"
            return 1
            
    # Construct the block of correlators.
    correlators = []
    for f in files:
        print f
        correlators.append(pion_correlator(f))
    correlators = np.array(correlators)
    
    # Basic checks on the output.
    print correlators.shape
    assert (len(files), nt) == correlators.shape
    
    # Write output.
    np.save(correlator_name(mspec), correlators)
    
    return 0
         
if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

