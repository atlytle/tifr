"""Delta_mix determination using \delta(m_v) = m^2_{vs} -m^2_{ss}/2."""

import sys
import numpy as np
import pylab as p

from HISQ_pion import pion as HHpion
from analyze_mixed_pions import pion as HOpion
from correlators import fit_cfuns, fit_single_cosh_osc
from correlators import fit_twopoint as fit_single_cosh
from soton.fits import line_fit2

sroot = '/Users/atlytle/Dropbox/TIFR/delta_mix/figs/'

def plot_dmsq(HOpions, HHpions, title=None, save=False, name=''):
    "Plot m^2_{vs} - m^2_{ss}/2."
    
    # Set up figure.
    fig = p.figure()
    p.rc('text', usetex=True)
    p.rc('font', size=16)
    p.rc('axes', linewidth=0.5)
    p.xlabel('$am_{v}$')
    p.ylabel('$m^2_{vs} - m^2_{ss}/2$')
    legend = ()

    xr = np.linspace(0.0,0.079)

    # First data set.
    hopions = HOpions[:4]
    r = HHpions[0]
    xs = [q.m2 for q in hopions]
    ys = [(q.msq - r.msq/2) for q in hopions]
    es = [nerror(q.sig_msq, r.sig_msq) for q in hopions]
    fit = line_fit2(zip(xs, ys, es))
    legend += p.errorbar(xs, ys, es, fmt='o', mfc='none', mec='b', lw=1)[0],
    # Fit results.
    p.errorbar(xr, fit.a + fit.b*xr, fmt='-b')
    plot = p.errorbar([0], [fit.a], [fit.sig_a], fmt='ob', lw=1, clip_on=False)
    # Adjust error bar.
    for x in plot[1]:
        x.set_clip_on(False)
        x.set_markeredgewidth(2)

    # Second data set.
    hopions = HOpions[4:]
    r = HHpions[1]
    xs = [q.m2 for q in hopions]
    ys = [(q.msq - r.msq/2) for q in hopions]
    es = [nerror(q.sig_msq, r.sig_msq) for q in hopions]
    fit = line_fit2(zip(xs, ys, es))
    legend += p.errorbar(xs, ys, es, fmt='s', mfc='none', mec='g', lw=1)[0],
    # Fit results.
    p.errorbar(xr, fit.a + fit.b*xr, fmt='-g')
    plot = p.errorbar([0], [fit.a], [fit.sig_a],
                      fmt='sg', lw=1, clip_on=False, alpha=0.7)
    # Adjust error bar.
    for x in plot[1]:
        x.set_clip_on(False)
        x.set_markeredgewidth(2)

    lg = p.legend(legend, ('$am_s = 0.0102$', '$am_s=0.0509$'), 'best', numpoints=1)
    lg.get_frame().set_linewidth(0.5)
    if save:
        p.savefig(name)
    else:
        p.show()

def nerror(sig1, sig2):
  "Naive error on delta m^2"
  return np.sqrt(sig1**2 + (sig2/2)**2)


def main(argv):
    # Load correlators.
    HHpions = [HHpion('0102', '0102'), 
               HHpion('0509', '0509')]

    HOpions = [HOpion('0102', '0165'),
               HOpion('0102', '0240'),
               HOpion('0102', '0380'), 
               HOpion('0102', '0731'), 
               HOpion('0509', '0165'),
               HOpion('0509', '0240'),
               HOpion('0509', '0380'), 
               HOpion('0509', '0731')] 

    # Calculate pion masses.
    for p in HHpions:
        print p.m1, p.m2
        f = fit_cfuns(p.JKcorrelators.real/100000, 10, 22, 64, fit_single_cosh)
        p.fit = f
        #print p.fit
        p.msq = (p.fit[0][1])**2
        p.sig_msq = 2*p.fit[1][1] # Naive error.
    
    for p in HOpions:
        print p.m1, p.m2
        f = fit_cfuns(p.JKcorrelators.real, 15, 29, 64, fit_single_cosh_osc)
        p.fit = f
        #print p.fit
        p.msq = (p.fit[0][2])**2
        p.sig_msq = 2*p.fit[1][2]  # Naive error.

    # Results.
    print ''
    r = HHpions[0]
    for p in HOpions[:4]:
        print p.m2, p.msq - r.msq/2, nerror(p.sig_msq, r.sig_msq)

    r = HHpions[1]
    for p in HOpions[4:]:
        print p.m2, p.msq - r.msq/2, nerror(p.sig_msq, r.sig_msq)

    # Plot results.
    plot_dmsq(HOpions, HHpions, save=True, name=sroot+'delta_msq.pdf')

    
    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
