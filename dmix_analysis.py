"""Delta_mix determination using \delta(m_v) = m^2_{vs} -m^2_{ss}/2."""

import sys
import numpy as np
import pylab as p

from HISQ_pion import pion as HHpion
from analyze_mixed_pions import pion as HOpion
from correlators import fit_cfuns, fit_single_cosh_osc
from correlators import fit_twopoint as fit_single_cosh


def main(argv):
    HHpions = [HHpion('0102', '0102'), 
               HHpion('0509', '0509')]

    HOpions = [HOpion('0102', '0165'), 
               HOpion('0102', '0380'), 
               HOpion('0102', '0731'), 
               HOpion('0509', '0165'), 
               HOpion('0509', '0380'), 
               HOpion('0509', '0731')] 

    for p in HHpions:
        print p.m1, p.m2
        f = fit_cfuns(p.JKcorrelators.real/100000, 10, 22, 64, fit_single_cosh)
        p.fit = f
        print p.fit
        p.msq = (p.fit[0][1])**2
    
    for p in HOpions:
        print p.m1, p.m2
        f = fit_cfuns(p.JKcorrelators.real, 15, 29, 64, fit_single_cosh_osc)
        p.fit = f
        print p.fit
        p.msq = (p.fit[0][2])**2

    print ''
    for p in HOpions[:3]:
        print p.m2, p.msq - HHpions[0].msq/2

    for p in HOpions[3:]:
        print p.m2, p.msq - HHpions[1].msq/2

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))