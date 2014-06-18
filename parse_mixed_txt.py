"""Parse pion data out of the correlator text files from Subhasis.

Returns the pion data as text files, as opposed to .npy format as in
in parse_mixed.py.  This makes it easier to do things like correlated
jackknife fits.
"""

import sys
import re

nt=64

ovpi_offset = 2626

def in_name(m1, m2, config):
    """Input text file."""
    root = '/Users/atlytle/Documents/dmix/for_tifr/corr/'
    return root + "mix_ksm{0}_ovm{1}_corr.{2}".format(m1, m2, config)

def out_name(m1, m2, config, type):
    root = '/Users/atlytle/Documents/dmix/pions/'
    if type == 'mix':
        return root + type + "_pion_m{0}_m{1}.{2}".format(m1, m2, config)
    elif type == 'hisq':
        return root + type + "_pion_m{0}.{1}".format(m1, config)
    elif type == 'ov':
        return root + type + "_pion_m{0}.{1}".format(m2, config)
    else:
        raise Exception('Pion type not recognized.')


def parse_pion(m1, m2, config, type):
    """Parse pion data. 

    Type = mix, hisq, or ov.
    """
    
    # Where to start reading data.
    if type == 'mix':
        offset = 28
    elif type == 'hisq':
        offset = 1327
    elif type == 'ov':
        offset = 2626
    else:
        raise Exception('Pion type not recognized.')
    
    # Extract pion data.
    result = []
    f = open(in_name(m1, m2, config), 'r')
    lines = f.readlines()
    for line in lines[offset:offset+nt]:
        result.append(line.rstrip())
    f.close()

    # Write the pion correlator.
    with open(out_name(m1, m2, config, type), 'w') as f:
        for line in result:
            f.write(line + '\n')

    return result

def main(argv):
    # Extract params from filename.
    s = re.compile('.+mix_ksm(\d+)_ovm(\d+)_corr\.(\d+)')
    
    # Validate input glob.
    for arg in argv:
        m = s.match(arg)
        assert m
        #m1, m2, config = m.group(1), m.group(2), m.group(3)
        #print m1, m2, config
    for arg in argv:
        m = s.match(arg)
        m1, m2, config = m.group(1), m.group(2), m.group(3)
        print m1, m2, config
        parse_pion(m1, m2, config, type='mix')
        parse_pion(m1, m2, config, type='hisq')
        parse_pion(m1, m2, config, type='ov')



    #print in_name('0509', '0165', '1000')
    #parse_pion('0509', '0165', '1000', type='mix')
    #parse_pion('0509', '0165', '1000', type='hisq')
    #parse_pion('0509', '0165', '1000', type='ov')



if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))