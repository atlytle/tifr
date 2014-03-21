import sys
import numpy as np

from read_overlap import correlator_name, correlator_name2
from resample import yield_JK_sample_np, JK_block, block
from calc_overlap import plot_correlator, plot_effmass
from correlators import fit_twopoint_cfuns

class pion:
    "Wrapper for pion correlation functions."
    def __init__(self, m1, m2):
        self.m1, self.m2 = m1, m2
        hbarc = 0.197327  # [GeV fm]
        self.ainv = hbarc/0.122  # [GeV], cf p.16 of 1212.4768
        self.correlators = np.load(correlator_name2(m1, m2))
        self.JKcorrelators = JK_block(self.correlators)
        
def main(argv=None):
    root = ''
    mlist = [0.0745, 0.001] #0.001 has no signal

    cs = [np.load(correlator_name(m)) for m in mlist]
    csJK = map(JK_block, cs)
    
    pions = [pion(0.001, 0.001), pion(0.0745, 0.0745)]
    
    #print cs[0].shape
    #print csJK[0].shape
    t=0
    for x in cs[0][0]:  # correlator
        print t, x
        t += 1
        
    print ''
    
    t=0
    for x in pions[1].correlators[0]:  # correlator
        print t, x
        t += 1

    #plot_correlator(csJK[1])
    #plot_effmass(csJK[1])
    #print fit_twopoint_cfuns(csJK[0]/10000., 10, 22, 64)

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))


