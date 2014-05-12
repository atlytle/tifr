"""Analysis of mixed pion correlators."""

import sys 
import numpy as np
import pylab as p

from parse_mixed import out_name
from resample import JK_block, JKsigma
#from calc_overlap import plot_correlator

class pion:
    "Wrapper for pion correlation functions."
    def __init__(self, m1, m2):
        self.m1s, self.m2s = m1, m2
        self.m1, self.m2 = float('0.'+m1), float('0.'+m2)  # Convert to float.
        hbarc = 0.197327  # [GeV fm]
        self.ainv = hbarc/0.122  # [GeV], cf p.16 of 1212.4768
        self.correlators = np.load(out_name(m1, m2))
        self.JKcorrelators = JK_block(self.correlators)

def plot_correlator(cfnc, save=False, name='', title=None):    
    ave, sigma = cfnc[0], JKsigma(cfnc)
    p.figure()
    if title is not None:
        p.title(title)
    p.xlabel('$t$')
    p.ylabel('$C[t]$')
    p.errorbar(range(len(ave)), ave.real,sigma.real, fmt='k-')
    if save:
        p.savefig(name)
    else:
        p.show()

def main(argv=None):
    pions = [pion('0102', '0731'), pion('0509', '0731'), pion('635', '0731')]
    
    root = '/Users/atlytle/Dropbox/TIFR/delta_mix/figs/'
    for p in pions:
        tag = 'HOpion_m{0}_m{1}'.format(p.m1s, p.m2s)
        sname = root + tag + '.pdf'
        plot_correlator(-p.JKcorrelators.real, True, sname, title=tag)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))