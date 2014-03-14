import sys
import struct
import numpy as np

from os.path import getsize

nx, ny, nz, nt = 24, 24, 24, 64
nc = 3
nfloat = nx*ny*nz*nt*nc*nc*2  # Number of 4byte numbers expected.

def correlator_name(m):
    return 'HHpion_l2464_m{0}_m{0}.npy'.format(m)

def correlator_name2(m1, m2):
    return 'HHpion_l2464_m{0}_m{1}.npy'.format(m1, m2)
    
def propagator_name(m, config):
    root = '/user2/atlytle/staggered'
    return root + '/hisq_{0}/prop_bin_hisq_{0}.{1}'.format(m, config)
    
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
    
def convert_to_complex(data):
    "Convert values to complex numbers."
    assert len(data)%2 == 0
    re, im = data[::2], data[1::2]
    return np.array(map(complex, re, im))
    
def pion_correlator2(file1, file2):
    "Construct pion correlator from two propagators."
    
    # Check file size.
    if getsize(file1) != (nfloat+31)*4:
        raise Exception('{0} does not have the expected size.'.format(file1))
    if getsize(file2) != (nfloat+31)*4:
        raise Exception('{0} does not have the expected size.'.format(file2))
    
    correlator = np.zeros((nt), dtype=complex)
    
    f1, f2 = open(file1, 'rb'), open(file2, 'rb')
    
    for f in f1, f2:
        # Read headers.  This consists of 22 4byte values.
        data = np.fromfile(f, dtype='<i', count=5)  # 5 floats.
        timestamp = np.fromfile(f, dtype='<c', count=64)  # 16 floats.
        order = np.fromfile(f, dtype='<i', count=1)  # 1 float.
    
    result = np.zeros((nt))
    
    for c in range(nc):
        color1 = np.fromfile(f1, dtype='<i', count=1)
        some_nums1 = np.fromfile(f1, dtype='<f', count=2) 
        data1 = np.fromfile(f1, dtype='<f', count=nfloat/3)
        color2 = np.fromfile(f2, dtype='<i', count=1)
        some_nums2 = np.fromfile(f2, dtype='<f', count=2) 
        data2 = np.fromfile(f2, dtype='<f', count=nfloat/3)
        for t in range(nt):
            c1 = convert_to_complex(extract_t(data1, t))
            c2 = convert_to_complex(extract_t(data2, t))
            correlator[t] += (c1*np.conj(c2)).astype(np.complex128).sum()
            
    f1.close()
    f2.close()
    
    return correlator
    
def convert_single_propagators(files):
    "Construct pion correlators from individual propagators."
    
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
            
                
def main():

    # Specify propagators to parse.
    config_list = \
    [1000, 1020, 1040, 1100, 1120, 1140, 1200, 1220, 1240, 1300, 1320, 1340,
     1400, 1420, 1440, 1500, 1520, 1540, 1600, 1620, 1640, 1700, 1720, 1740,
     1800, 1820, 1840, 1900, 1920, 1940, 2000, 2020, 2040, 2100, 2120, 2140,
     2200, 2220, 2240, 2300, 2320, 2340, 2400, 2420, 2440, 2500, 2520, 2540,
     2600, 2620, 2640, 2700, 2720, 2740, 2800, 2820, 2840, 2900, 2920, 2940,
     3000, 3020, 3040, 3100, 3120, 3140, 3200, 3220, 3240, 3300, 3320, 3340,
     3400, 3420, 3440, 3500, 3600, 3700, 3720, 3800, 3900]
    
    flist1 = [propagator_name('635', config) for config in config_list]
    flist2 = [propagator_name('0509', config) for config in config_list]
    
    # Construct the block of correlators.
    correlators = []
    for f1, f2 in zip(flist1, flist2):
        print f1
        print f2, '\n'
        correlators.append(pion_correlator2(f1, f2))
    correlators = np.array(correlators)
    
    # Basic checks on the output.
    print correlators.shape
    assert (len(flist1), nt) == correlators.shape
    
    # Write output.
    np.save(correlator_name2('635', '0509'), correlators)
  
    return 0
         
if __name__ == "__main__":
    sys.exit(main())

