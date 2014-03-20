import sys
import numpy as np
import pylab as p

from read_HISQ import correlator_name, correlator_name2
from resample import JK_block, JKsigma
from calc_overlap import plot_correlator #fit_twopoint_cfuns,plot_effmass
from correlators import fit_twopoint_cfuns, fit_cfuns_double_cosh_osc, naive_effmass

def plot_effmass(cfnc, save=False, name='', title=None, yrange=[]):
    effmass = naive_effmass(cfnc)
    ave, sigma = effmass[0], JKsigma(effmass)
    p.figure()
    if title is not None:
        p.title(title)
    if yrange:
        p.ylim(yrange)
    p.xlabel('$t$')
    p.ylabel('$|\\frac{C[t+1]}{C[t]}|$')
    p.errorbar(range(len(ave)), ave, sigma, fmt='ko')
    if save:
        p.savefig(name)
    else:
        p.show()
        
def plot_effmass_with_fit(cfnc, ti, tf, T, 
                          save=False, name='', title=None, yrange=[]):
    "Also show the mass determined from a fit."
    effmass = naive_effmass(cfnc)
    ave, sigma = effmass[0], JKsigma(effmass)
    fit = fit_twopoint_cfuns(cfnc, ti, tf, T)
    m, sigm = fit[0][1], fit[1][1]
    p.figure()
    if title is not None:
        p.title(title)
    if yrange:
        p.ylim(yrange)
    p.xlabel('$t$')
    p.ylabel('$|\\frac{C[t+1]}{C[t]}|$')
    p.errorbar(range(len(ave)), ave, sigma, fmt='ko')
    p.hlines([m-sigm, m+sigm], ti, tf, colors='r')
    #p.axhspan((m-sigm), (m+sigm), facecolor=(1.0,0.0,0.0), alpha=0.5) 
    if save:
        p.savefig(name)
    else:
        p.show()
        
class pion:
    "Wrapper for pion correlation functions."
    def __init__(self, m1, m2):
        self.m1, self.m2 = m1, m2
        self.correlators = np.load(correlator_name2(m1, m2))
        self.JKcorrelators = JK_block(self.correlators)

def main(argv=None):
    #root = ''
    #mlist = ['635', '0102', '0509']
    #cs = [np.load(correlator_name(m)) for m in mlist]
    #csJK = map(JK_block, cs)

    #c_635_0509 = np.load(correlator_name2('635', '0509'))
    #c_635_0509JK = JK_block(c_635_0509)


    c_0509_635 = np.load(correlator_name2('0509', '635'))
    c_0509_635JK = JK_block(c_0509_635)

    c_0509_0509 = np.load(correlator_name2('0509', '0509'))
    c_0509_0509JK = JK_block(c_0509_0509)
    
    test = pion('0509', '0509')
    
    assert (test.JKcorrelators == c_0509_0509JK).all()


    print c_0509_0509.shape

    #t = 0
    #for x in c_0509_635JK[0]:
    #    print t, x.real,' ', x.imag
    #    t += 1
        
        
    plot_correlator(c_0509_635JK.real)
    plot_effmass(c_0509_635JK.real)
    print fit_twopoint_cfuns(c_0509_635JK.real/10000., 10, 22, 64)
    print fit_cfuns_double_cosh_osc(c_0509_635JK.real/10000., 10, 22, 64)
    plot_effmass_with_fit(c_0509_635JK.real, 10, 22, 64)

#    for ti in range(8, 15):
#        for tf in range(18, 32):
#            print ti, tf, ": ", fit_twopoint_cfuns(c_0509_635JK.real/10000., ti, tf, 64)



    #plot_correlator(c_0509_635JK.real)
    #plot_effmass(c_0509_635JK.real)
    #print fit_twopoint_cfuns(c_0509_635JK.real/10000., 10, 22, 64)

    #plot_correlator(csJK[0])
    #plot_correlator(csJK[1])
    #plot_correlator(csJK[2])

    #plot_effmass(csJK[0])
    #plot_effmass(csJK[1])
    #plot_effmass(csJK[2])

    #print fit_twopoint_cfuns(csJK[0]/10000., 10, 22, 64)
    #print fit_twopoint_cfuns(csJK[1]/100000., 10, 22, 64)
    #print fit_twopoint_cfuns(csJK[2]/10000., 10, 22, 64)

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
