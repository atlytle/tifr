"""Analysis of mixed pion correlators."""

import sys 
import numpy as np
import pylab as p
from math import cosh

from parse_mixed import pion_name, rhox_name, rho_name, scalar_name
from resample import JK_block, JKsigma
from correlators import fit_cfuns_double_cosh_osc, fit_cfuns, fit_single_cosh_osc
#from calc_overlap import plot_correlator


class corr:
    def __init__(self, m1, m2):
        self.m1s, self.m2s = m1, m2
        self.m1, self.m2 = float('0.'+m1), float('0.'+m2)  # Convert to float.
        hbarc = 0.197327  # [GeV fm]
        self.ainv = hbarc/0.122  # [GeV], cf p.16 of 1212.4768
        self.correlators = np.load(self.name)
        self.JKcorrelators = JK_block(self.correlators)

class pion(corr):
    "Wrapper for pion correlation functions."
    def __init__(self, m1, m2):
        self.name = pion_name(m1, m2)
        corr.__init__(self, m1, m2)

class rhox(corr):
    "Wrapper for rhox correlation functions."
    def __init__(self, m1, m2):
        self.name = rhox_name(m1, m2)
        corr.__init__(self, m1, m2)

class rho(corr):
    "Wrapper for rho correlation functions."
    def __init__(self, m1, m2):
        self.name = rho_name(m1, m2)
        corr.__init__(self, m1, m2)

class scalar(corr):
    "Wrapper for scalar correlation functions."
    def __init__(self, m1, m2):
        self.name = scalar_name(m1, m2)
        corr.__init__(self, m1, m2)

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

def plot_masses(ps, rs, scs, save=False, name=''):
    p.figure
    p.rc('text', usetex=True)
    p.rc('font', size=18)
    p.xlabel('$am_1 + am_2$')
    p.ylabel('$aM_{vs}$')
    legend = ()
    
    # Pions.
    xs, ys, es = zip(*[q.pt for q in ps])  # Unpack data.
    legend += p.errorbar(xs, ys, es, fmt='o')[0],

    # Rhoxs.
    #xs, ys, es = zip(*[r.pt for r in rxs])  # Unpack data.
    #legend += p.errorbar(xs, ys, es, fmt='o')[0],

    # Rhos.
    xs, ys, es = zip(*[r.pt for r in rs])  # Unpack data.
    legend += p.errorbar(xs, ys, es, fmt='o')[0],

    # Scalars.
    xs, ys, es = zip(*[r.pt for r in scs])  # Unpack data.
    legend += p.errorbar(xs, ys, es, fmt='o')[0],
    
    p.legend(legend, ('$\pi$', r'$\rho$', '$a_0$'), 'best')
    if save:
        p.savefig(name)
    else:
        p.show()

fitfunc = lambda p, x: p[0]*cosh(p[1]*(x-32.)) + \
                          ((-1)**x)*p[2]*cosh(p[3]*(x-32.))

def main(argv=None):
    pions = [pion('0102', '0731'), pion('0509', '0731'), 
             pion('0102', '0165'), pion('0509', '0165')]
    rhoxs = [rhox('0102', '0731'), rhox('0509', '0731'),
             rhox('0102', '0165'), rhox('0509', '0165')]
    rhos = [rho('0102', '0731'), rho('0509', '0731'),
            rho('0102', '0165'), rho('0509', '0165')]
    scalars = [scalar('0102', '0731'), scalar('0509', '0731'),
               scalar('0102', '0165'), scalar('0509', '0165')]


    fit = fit_cfuns_double_cosh_osc
    
    
    root = '/Users/atlytle/Dropbox/TIFR/delta_mix/figs/'
    '''
    print "pions"
    for p in pions:
        #tag = 'HOpion_m{0}_m{1}'.format(p.m1s, p.m2s)
        #sname = root + tag + '.pdf'
        #plot_correlator(-p.JKcorrelators.real, False, sname, title=tag)
        print "double"
        print fit(-p.JKcorrelators.real, 14, 29, 64)
        print "single"
        print fit_cfuns(-p.JKcorrelators.real, 14, 29, 64, fit_single_cosh_osc)
        #print -p.JKcorrelators.real[0
        print ''
    #print [fitfunc(v[0], x) for x in range(64)]

    print "rhos"
    for r in rhos:
        print "double"
        print fit(-r.JKcorrelators.real, 14, 29, 64)
        print "single"
        print fit_cfuns(-r.JKcorrelators.real, 14, 29, 64, fit_single_cosh_osc)
        print ''
    '''

    for p in pions:
        f = fit_cfuns(p.JKcorrelators.real, 15, 29, 64, fit_single_cosh_osc)
        p.pt = (p.m1+p.m2, f[0][2], f[1][2])
        #print p.pt

    for r in rhoxs:
        f = fit_cfuns(r.JKcorrelators.real, 15, 29, 64, fit_single_cosh_osc)
        r.pt = (r.m1+r.m2, f[0][2], f[1][2])
        #print r.pt

    for r in rhos:
        f = fit_cfuns(r.JKcorrelators.real, 15, 29, 64, fit_single_cosh_osc)
        r.pt = (r.m1+r.m2, f[0][2], f[1][2])
        #print r.pt

    for r in scalars:
        f = fit_cfuns(r.JKcorrelators.real, 15, 29, 64, fit_single_cosh_osc)
        r.pt = (r.m1+r.m2, f[0][2], f[1][2])
        #print r.pt


    sname = root +'mixed_mpi2.pdf'
    plot_masses(pions, rhos, scalars, save=False, name=sname)
    #print rhos[0].JKcorrelators.real[0]
    #print pions[0].JKcorrelators.real[0]
    #print fit(rhos[0].JKcorrelators.real, 10, 29, 64)
    #print fit_cfuns(rhos[0].JKcorrelators.real, 10, 29, 64, fit_single_cosh_osc)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))