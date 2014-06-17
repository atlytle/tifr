import sys
import numpy as np

from read_overlap import correlator_name, correlator_name2
from resample import yield_JK_sample_np, JK_block, block
from calc_overlap import plot_correlator, plot_effmass
from correlators import fit_cfuns, fit_twopoint as fit_single_cosh

class pion:
    "Wrapper for pion correlation functions."
    def __init__(self, m1, m2):
        self.m1, self.m2 = m1, m2
        hbarc = 0.197327  # [GeV fm]
        self.ainv = hbarc/0.122  # [GeV], cf p.16 of 1212.4768
        root = '/Users/atlytle/Dropbox/pycode/tifr/data/'
        self.correlators = np.load(root+correlator_name2(m1, m2))
        self.JKcorrelators = JK_block(self.correlators)

def plot_Msq_over_m(dat, save=False, name=''):
    pass

        
def main(argv=None):
    
    pions = [pion(0.001, 0.001), pion(0.009, 0.009), pion(0.0165, 0.0165),
             pion(0.024, 0.024), pion(0.038, 0.038)]
    
    
    # m=0.001 seems unstable    
    for p in pions[1:]:
        print p.m1, p.m2
        f = fit_cfuns(p.JKcorrelators.real/10000, 8, 24, 64, fit_single_cosh)
        #print f
        print f[0][1], f[1][1], f[2]
        print ''


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))


