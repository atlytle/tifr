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
        
def plot_Msq_over_msum(dat, save=False, name=''):
    "Plot (msum, Msq/msum, sig) data."
    msum, Msq_over_msum, sig = zip(*dat)
    p.figure()
    p.title('HH pseudoscalars - Preliminary')
    p.xlabel("$m_1 + m_2$")
    p.ylabel('$(a M_\pi)^2/a(m_1 + m_2)$')
    p.errorbar(msum, Msq_over_msum, sig, fmt='ko')
    if save:
        p.savefig(name)
    else:
        p.show()
    
        
class pion:
    "Wrapper for pion correlation functions."
    def __init__(self, m1, m2):
        self.root = '/Users/atlytle/Dropbox/pycode/tifr/data/'
        self.m1, self.m2 = float('0.'+m1), float('0.'+m2)  # Convert to float.
        hbarc = 0.197327  # [GeV fm]
        self.ainv = hbarc/0.122  # [GeV], cf p.16 of 1212.4768
        self.correlators = np.load(self.root + correlator_name2(m1, m2))
        self.JKcorrelators = JK_block(self.correlators)

def main(argv=None):

    root = '/Users/atlytle/Dropbox/TIFR/delta_mix/figs/'

    
    pions = [pion('0102', '0102'), pion('0102', '0509'), pion('0509', '0509'),
             pion('0102', '635'), pion('0509', '635'), pion('635', '635')]
             
    for p in pions:
        print p.m1, p.m2
        p.msum = p.m1 + p.m2
        fit = fit_twopoint_cfuns(p.JKcorrelators.real/100000, 10, 22, 64)
        p.M, p.sigM, p.chi2 = fit[0][1], fit[1][1], fit[2]
        print p.msum*p.ainv, (p.M)*(p.ainv)
        print ''
        
    mass_dat = [(p.msum, (p.M**2)/(p.msum), 2*p.sigM/p.msum) for p in pions]  # Naive error.
    plot_Msq_over_msum(mass_dat, save=True, name=root+'HH_pseudo.pdf')
        
    # Getting a horrible chi^2:
    print pions[3].m1, pions[3].m2
    plot_effmass(pions[3].JKcorrelators.real)
        

        
#    plot_correlator(c_0509_635JK.real)
#    plot_effmass(c_0509_635JK.real)
#    print fit_twopoint_cfuns(c_0509_635JK.real/10000., 10, 22, 64)
#    print fit_cfuns_double_cosh_osc(c_0509_635JK.real/10000., 10, 22, 64)
#    plot_effmass_with_fit(c_0509_635JK.real, 10, 22, 64)

#    for ti in range(8, 15):
#        for tf in range(18, 32):
#            print ti, tf, ": ", fit_twopoint_cfuns(c_0509_635JK.real/10000., ti, tf, 64)

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
