import sys
import struct
import itertools
import numpy as np

from os.path import getsize

nx, ny, nz, nt = 24, 24, 24, 64
V = nx*ny*nz  # Spatial volume.
nc = 3  # N_c.
ns = 4  # N_spin.
ndouble = V*nt*nc*nc*ns*ns*2  # Number of 8byte numbers expected.

def prop_name(config, m):
    root = '/user5/nilmani/projects/HISQ/quark_prop/L24T64/'
    return root + 'l2464f211b600m0102m0509m635a.'\
    '{0}.gf_hyp_prop_src_wall_t00_m{1:7f}'.format(config, m)

def correlator_name(m):
    return 'OOpion_l2464__m{0:7f}_m{0:7f}.npy'.format(m)
    
def correlator_name2(m1, m2):
    return 'OOpion_l2464__m{0:7f}_m{1:7f}.npy'.format(m1, m2)

def pion_t(filename, t):
    '''Extract the pion correlation function at time t from the propagator.'''
    tmp = []
    for i in range(2*nc*nc*ns*ns):
        with open(filename, "rb") as f:  # Inefficient?
            f.seek(i*8*nt*V + 8*t*V, 0)
            data = np.fromfile(f, dtype='>d', count=V)
            tmp.append(np.sum(np.square(data)))
    return np.sum(tmp)
    
def extract_t5(filename, t):
    '''Extract data at timeslice t from the propagator.

    Converts the raw data into complex numbers.
    '''
    # Loop structure: s c r/i s c t z y x.
    # Pluck out bits at t. Store in tmp.
    tmp = np.array([])
    for i in range(2*nc*nc*ns*ns):
        with open(filename, "rb") as f:  # Inefficient?
            f.seek(i*8*nt*V + 8*t*V, 0)
            data = np.fromfile(f, dtype='>d', count=V)
            tmp = np.hstack((tmp, data))
    # Convert to complex numbers.  Store in ctmp.
    ctmp = np.array([])
    for chunk in np.hsplit(tmp, ns*nc):
        chunk_re, chunk_im = np.hsplit(chunk, 2)
        chunk_c = map(complex, chunk_re, chunk_im)
        ctmp = np.hstack((ctmp, chunk_c))
        
    return ctmp
            
def pion_correlator(filename):
    "Construct pion correlator from propagator."
    correlator = np.zeros((nt))
    for t in range(nt):
        correlator[t] = pion_t(filename, t)
    return correlator
    
def pion_correlator2(file1, file2):
    "Construct pion correlator from two propagators."
    correlator = np.zeros((nt), dtype=complex)
    for t in range(nt):
        nums1 = extract_t5(file1, t)
        nums2 = extract_t5(file2, t)
        correlator[t] = np.sum(nums1*np.conj(nums2))
    return correlator

def check_length(filename):
    '''Ensure the file has the expected size.'''
    if getsize(filename) != ndouble*8:
        raise Exception('{0} does not have the expected size.'.format(filename))

def convert_single_propagators(files):
    "Construct pion correlators from individual overlap propagators."
    # Some basic checks on the input.
    head0, config0, middle0, mass0 = files[0].split('.')
    for f in files:
        check_length(f)
        head, config, middle, mass = f.split('.')
        if (head != head0) or (middle != middle0) or (mass != mass0):
            print "You might not want to combine these!"
            return 1

    # Construct the block of correlators.
    correlators = []
    for f in files:
        correlators.append(pion_correlator(f))
    correlators = np.array(correlators)

    # Basic checks on the output.
    print correlators.shape
    assert (len(files), nt) == correlators.shape

    # Write output.
    m = float('0.'+mass0)
    np.save(correlator_name(m), correlators) 

def main(files):
    f1, f2 = files[0], files[1]
    
    x = pion_correlator2(f1, f2)
    t = 0
    for c in x:
        print t, c.real, c.imag
        t += 1
    
    return 0
    
#    # Specify propagators to parse.
#    config_list = [1020, 1040]
#    flist1 = [propagator_name(config, 0.0745) for config in config_list]
#    flist2 = [propagator_name(config, 0.001) for config in config_list]
#    
#    # Construct the block of correlators.
#    correlators = []
#    for f1, f2 in zip(flist1, flist2):
#        print f1
#        print f2, '\n'

    
    return 0
   
         
if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

