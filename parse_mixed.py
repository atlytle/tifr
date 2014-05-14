"""Parse mixed correlator output files into numpy arrays."""

import re
import sys
import numpy as np

def in_name(m1, m2, config):
    """Input text file."""
    root = "/Users/atlytle/Dropbox/pycode/tifr/data/for_tifr/corr/"
    return root + "mix_ksm{0}_ovm{1}_corr.{2}".format(m1, m2, config)

def out_name(m1, m2):
    """Output numpy file."""
    root = "/Users/atlytle/Dropbox/pycode/tifr/data/"
    return root + "HOpion_l2464_m{0}_m{1}.npy".format(m1, m2)
pion_name = out_name

def rhox_name(m1, m2):
    """Output rho_x numpy file."""
    root = "/Users/atlytle/Dropbox/pycode/tifr/data/"
    return root + "HOrhox_l2464_m{0}_m{1}.npy".format(m1, m2)

def rho_name(m1, m2):
    """Output rho numpy file."""
    root = "/Users/atlytle/Dropbox/pycode/tifr/data/"
    return root + "HOrho_l2464_m{0}_m{1}.npy".format(m1, m2)

def scalar_name(m1, m2):
    """Output scalar numpy file."""
    root = "/Users/atlytle/Dropbox/pycode/tifr/data/"
    return root + "HOscalar_l2464_m{0}_m{1}.npy".format(m1, m2)

def parse_pion(m1, m2, config):
    """Parse the pion data and convert to numpy array."""
    f = open(in_name(m1, m2, config), 'r')
    
    # Skip to the mixed pion correlator section.
    # This assumes it is the first PION entry.
    x = f.readline()
    while x:
        if re.match('correlator:\s+PION', x):
            break
        x = f.readline()

    # Throw away header.
    print x
    for i in range(5):
        print f.readline().strip()

    result = []
    for i in range(64):
        t, r, im = f.readline().strip().split('\t')
        result.append(complex(float(r), float(im)))        
    
    f.close()

    return np.array(result)

def parse_rhox(m1, m2, config):
    """Parse the rho_x data and convert to numpy array."""
    f = open(in_name(m1, m2, config), 'r')
    
    # Skip to the mixed rho correlator section.
    # This assumes it is the first RHOX entry.
    x = f.readline()
    while x:
        if re.match('correlator:\s+RHOX', x):
            break
        x = f.readline()

    # Throw away header.
    print x
    for i in range(5):
        print f.readline().strip()

    result = []
    for i in range(64):
        t, r, im = f.readline().strip().split('\t')
        result.append(complex(float(r), float(im)))        
    
    f.close()

    return np.array(result)

def parse_rho(m1, m2, config):
    """Parse the rho data and convert to numpy array."""
    f = open(in_name(m1, m2, config), 'r')
    
    # Skip to the mixed rho correlator section.
    # This assumes it is the first RHO* entry. !!
    x = f.readline()
    while x:
        if re.match('correlator:\s+RHO', x):
            break
        x = f.readline()

    # Throw away header.
    print x
    for i in range(5):
        print f.readline().strip()

    result = []
    for i in range(64):
        t, r, im = f.readline().strip().split('\t')
        result.append(complex(float(r), float(im)))        
    
    f.close()

    return np.array(result)

def parse_scalar(m1, m2, config):
    """Parse the rho data and convert to numpy array."""
    f = open(in_name(m1, m2, config), 'r')
    
    # Skip to the mixed scalar correlator section.
    # This assumes it is the first SCALAR entry.
    x = f.readline()
    while x:
        if re.match('correlator:\s+SCALAR', x):
            break
        x = f.readline()

    # Throw away header.
    print x
    for i in range(5):
        print f.readline().strip()

    result = []
    for i in range(64):
        t, r, im = f.readline().strip().split('\t')
        result.append(complex(float(r), float(im)))        
    
    f.close()

    return np.array(result)

def main(argv):
    
    # Extract params from filename.
    s = re.compile('mix_ksm(\d+)_ovm(\d+)_corr\.(\d+)')
    
    # Validate input glob.
    m = s.match(argv[0])
    assert m
    m1, m2 = m.group(1), m.group(2)
    print m1, m2
    for arg in argv:
        m = s.match(arg)
        assert m
        assert m1 == m.group(1)
        assert m2 == m.group(2)
        
    # Construct numpy correlators.
    
    # Pion.
    result = []
    for arg in argv:
        m = s.match(arg)
        config = m.group(3)
        print config
        result.append(parse_pion(m1, m2, config))
    result = np.array(result)
    print result.shape
    np.save(pion_name(m1,m2), result)

    # Rho.
    result = []
    for arg in argv:
        m = s.match(arg)
        config = m.group(3)
        print config
        result.append(parse_rho(m1, m2, config))
    result = np.array(result)
    print result.shape
    np.save(rho_name(m1,m2), result)

    # Rho_x.
    result = []
    for arg in argv:
        m = s.match(arg)
        config = m.group(3)
        print config
        result.append(parse_rhox(m1, m2, config))
    result = np.array(result)
    print result.shape
    np.save(rhox_name(m1,m2), result)

    # Scalar.
    result = []
    for arg in argv:
        m = s.match(arg)
        config = m.group(3)
        print config
        result.append(parse_scalar(m1, m2, config))
    result = np.array(result)
    print result.shape
    np.save(scalar_name(m1,m2), result)

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

