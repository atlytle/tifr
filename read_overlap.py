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
    return 'l2464f211b600m0102m0509m635a.'\
    '{0}.gf_hyp_prop_src_wall_t00_m{1:7f}'.format(config, m)

def correlator_name(m):
    return 'OOpion_l2464__m{0:7f}_m{0:7f}.npy'.format(m)

def extract_t(filename, t):
    '''Extract the pion correlation function at time t from the propagator.'''
    tmp = []
    for i in range(2*nc*nc*ns*ns):
        with open(filename, "rb") as f:  # Inefficient.
            f.seek(i*8*nt*V + 8*t*V, 0)
            data = np.fromfile(f, dtype='>d', count=V)
            tmp.append(np.sum(np.square(data)))
    return np.sum(tmp)

def pion_correlator(filename):
    '''Construct pion correlator from propagator.'''
    correlator = np.zeros((nt))
    for t in range(nt):
        correlator[t] = extract_t(filename, t)
    return correlator

def check_length(filename):
    '''Ensure the file has the expected size.'''
    if getsize(filename) != ndouble*8:
        raise Exception('{0} does not have the expected size.'.format(filename))


def main(files):
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
    
    return 0
   
         
if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

