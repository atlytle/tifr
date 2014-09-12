import sys
import struct
import itertools
import numpy as np

from os.path import getsize

from gammas import *

ndouble = V*nt*nc*nc*ns*ns*2  # Number of 8byte numbers expected.

def prop_name(m, config):
    #root = '/user5/nilmani/projects/HISQ/quark_prop/L24T64/'
    root = '/Users/atlytle/Dropbox/pycode/tifr/check_Subhasis/'
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

ar = np.array
def extract_t5(filename, t):
    '''Extract data at timeslice t from the propagator.

    Converts the raw data into complex numbers.
    '''
    # Loop structure: s c r/i s c t z y x.
    # Pluck out bits at t. Store in tmp.
    tmp = []
    for i in range(2*nc*nc*ns*ns):
        with open(filename, "rb") as f:  # Inefficient?
            f.seek(i*8*nt*V + 8*t*V, 0)
            data = np.fromfile(f, dtype='>d', count=V)
            tmp.append(data)
    tmp = ar(tmp, dtype=np.float).reshape((-1,))
    # Convert to complex numbers.  Store in ctmp.
    ctmp_re = []
    ctmp_im = []
    for chunk in np.hsplit(tmp, ns*nc):
        chunk_re, chunk_im = np.hsplit(chunk, 2)
        ctmp_re.append(chunk_re)
        ctmp_im.append(chunk_im)
    ctmp_re = ar(ctmp_re, dtype=np.float).reshape((-1,))
    ctmp_im = ar(ctmp_im, dtype=np.float).reshape((-1,))
    ctmp = ctmp_re + 1j*ctmp_im
        
    return ctmp
            
def pion_correlator(filename):
    "Construct pion correlator from propagator."
    correlator = np.zeros((nt))
    for t in range(nt):
        print t
        correlator[t] = pion_t(filename, t)
    return correlator
    
def pion_correlator2(file1, file2):
    "Construct pion correlator from two propagators."
    correlator = np.zeros((nt), dtype=complex)
    for t in range(nt):
        print t
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
        print f
        correlators.append(pion_correlator(f))
    correlators = np.array(correlators)

    # Basic checks on the output.
    print correlators.shape
    assert (len(files), nt) == correlators.shape

    # Write output.
    m = float('0.'+mass0)
    np.save(correlator_name(m), correlators) 
    print correlators[0]

def main(files):
   
    #print pion_correlator2(files[0], files[0])
    convert_single_propagators(files)
    
    return 0
         
if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

