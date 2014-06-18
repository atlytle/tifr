"""Parse pion data out of the correlator text files from Subhasis.

Returns the pion data as text files, as opposed to .npy format as in
in parse_mixed.py.  This makes it easier to do things like correlated
jackknife fits.
"""

import sys
import re
import numpy as np

from resample import JK_block

nt=64

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
        raise Exception('Pion type {0} not recognized.'.format(type))


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

def read_pion(m1, m2, config_list, type):
    """Read pion data from individual txt files.

    Returns numpy correlator block ordered by config_list.
    Note that type hisq/ov only uses mass m1/m2.
    """
    result = []
    for config in config_list:
        corr = []
        f = open(out_name(m1, m2, config, type), 'r')
        for i in range(nt):
            t, re, im = f.readline().strip().split('\t')
            corr.append(complex(float(re),float(im)))
        result.append(corr)
        f.close()
    return np.array(result)

class corr:
    def __init__(self, m1, m2, config_list, type):
        self.type = type
        self.m1s, self.m2s = m1, m2
        self.m1, self.m2 = float('0.'+m1), float('0.'+m2)  # Convert to float.
        hbarc = 0.197327  # [GeV fm]
        self.ainv = hbarc/0.122  # [GeV], cf p.16 of 1212.4768
        self.correlators = read_pion(m1, m2, config_list, self.type)
        self.JKcorrelators = JK_block(self.correlators)

class HHpion(corr):
    def __init__(self, m1, config_list):
        self.type = 'hisq'
        corr.__init__(self, m1, m1, config_list, self.type)

class HOpion(corr):
    def __init__(self, m1, m2, config_list):
        self.type = 'mix'
        corr.__init__(self, m1, m2, config_list, self.type)

def main(argv):
    # # Extract params from filename.
    # s = re.compile('.+mix_ksm(\d+)_ovm(\d+)_corr\.(\d+)')
    
    # # Validate input glob.
    # for arg in argv:
    #     m = s.match(arg)
    #     assert m
    #     #m1, m2, config = m.group(1), m.group(2), m.group(3)
    #     #print m1, m2, config
    # for arg in argv:
    #     m = s.match(arg)
    #     m1, m2, config = m.group(1), m.group(2), m.group(3)
    #     print m1, m2, config
    #     parse_pion(m1, m2, config, type='mix')
    #     parse_pion(m1, m2, config, type='hisq')
    #     parse_pion(m1, m2, config, type='ov')
    
    config_list = [1000, 1020]
    p = read_pion('0102', '0380', config_list, 'mix')
    p2 = read_pion('0102', '0380', config_list, 'hisq')

    print p.shape
    print p[1]
    print p2[1]

    hh = HHpion('0102', config_list)
    ho = HOpion('0102', '0380', config_list)
    print hh.correlators[1]
    print ho.correlators[1]



if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))