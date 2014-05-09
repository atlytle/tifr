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
    result = []
    for arg in argv:
        m = s.match(arg)
        config = m.group(3)
        print config
        result.append(parse_pion(m1, m2, config))
    result = np.array(result)
    print result.shape
    np.save(out_name(m1,m2), result)

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

