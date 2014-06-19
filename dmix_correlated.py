"""Delta_mix determination using \delta(m_v) = m^2_{vs} -m^2_{ss}/2.

Correlated version of dmix_analysis.py.
"""

import sys
import numpy as np
import pylab as p

from parse_mixed_txt import HHpion, HOpion
from correlators import fit_cfuns, fit_cfuns_through, fit_single_cosh_osc
from correlators import fit_twopoint as fit_single_cosh
from resample import JKsigma
from soton.fits import line_fit2

sroot = '/Users/atlytle/Dropbox/TIFR/delta_mix/figs/'

list_0102_0165 = [1000, 1020, 1040, 1100, 1140, 1200, 1220, 1240, 1300, 1320,
                  1340, 1400, 1420, 1500, 1520, 1540, 1600, 1620, 1700, 1720,
                  1800, 1820, 1840, 1900, 1920, 2000, 2020, 2100, 2140, 2200,
                  2220, 2300, 2320, 2400, 2440, 2500, 2540, 2600, 2620, 2700,
                  3000, 3100, 3200, 3500]

list_0509_0165 = list_0102_0165  # Overlap configs were limiting factor.

list_0102_0240 = [1000, 1020, 1040, 1100, 1140, 1200, 1220, 1240, 1300, 1320,
                  1340, 1400, 1420, 1500, 1520, 1600, 1620, 1720, 2000, 2020,
                  2100, 2140, 2200, 2220, 2300, 2320, 2400, 2440, 2540, 2700,
                  2720, 2800, 2900, 3000, 3100, 3200, 3300, 3400, 3500, 3600,
                  3700, 3800, 3900]

list_0509_0240 = list_0102_0240

list_0102_0380 = [1000, 1020, 1040, 1100, 1140, 1200, 1220, 1240, 1300, 1320,
                  1340, 1400, 1420, 1500, 1520, 1540, 1600, 1620, 1800, 2000,
                  2020, 2100, 2140, 2200, 2220, 2300, 2320, 2400, 2440, 2500,
                  2540, 2620, 2700, 2800, 2900, 3100, 3200, 3600, 3900]

list_0509_0380 = list_0102_0380

list_0102_0731 = [1000, 1020, 1040, 1100, 1140, 1220, 1320, 1340, 2000, 2020,
                  2100, 2140, 2200, 2220, 2240, 2300, 2320, 2340, 2400, 2500,
                  2600, 2700, 2720, 2800, 2900, 3000, 3100, 3200, 3300, 3400,
                  3500, 3600, 3700, 3800, 3900]

list_0509_0731 = list_0102_0731

def dmsq(HOpion, HHpion):
    """(delta_m^2, sigma) from jackknife samples of m^2."""
    dmsqJK = HOpion.msqJK - HHpion.msqJK/2
    return dmsqJK[0], JKsigma(dmsqJK)

def plot_dmsq(HOpions, HHpions, title=None, save=False, name=''):
    "Plot m^2_{vs} - m^2_{ss}/2."
    
    # Set up figure.
    fig = p.figure()
    p.rc('text', usetex=True)
    p.rc('font', size=16)
    p.rc('axes', linewidth=0.5)
    p.xlabel('$m_{v} \, r_1$')
    p.ylabel('$\Biggl( \, m_{vs}^2 - m_{ss}^2/2 \, \Biggr) \cdot \, r_1^2$')
    legend = ()

    xr = np.linspace(0.0,0.2)
    r1oa = HOpions[0].r1oa  # Conversion factor r1/a.

    # First data set.
    hopions = HOpions[:4]
    hhpions = HHpions[:4]
    xs = np.array([q.m2 for q in hopions])*r1oa
    ys = np.array([dmsq(*x)[0] for x in zip(hopions, hhpions)])*r1oa*r1oa
    es = np.array([dmsq(*x)[1] for x in zip(hopions, hhpions)])*r1oa*r1oa
    print ys
    print es
    fit = line_fit2(zip(xs, ys, es))
    legend += p.errorbar(xs, ys, es, fmt='o', mfc='none', mec='b', lw=1)[0],
    # Fit results.
    p.errorbar(xr, fit.a + fit.b*xr, fmt='-b')
    plot = p.errorbar([0], [fit.a], [fit.sig_a], fmt='ob', lw=1, clip_on=False)
    print fit.a, fit.sig_a
    # Adjust error bar.
    for x in plot[1]:
        x.set_clip_on(False)
        x.set_markeredgewidth(2)

    # Second data set.
    hopions = HOpions[4:]
    hhpions = HHpions[4:]
    xs = np.array([q.m2 for q in hopions])*r1oa
    ys = np.array([dmsq(*x)[0] for x in zip(hopions, hhpions)])*r1oa*r1oa
    es = np.array([dmsq(*x)[1] for x in zip(hopions, hhpions)])*r1oa*r1oa
    print ys
    print es
    fit = line_fit2(zip(xs, ys, es))
    legend += p.errorbar(xs, ys, es, fmt='s', mfc='none', mec='g', lw=1)[0],
    # Fit results.
    p.errorbar(xr, fit.a + fit.b*xr, fmt='-g')
    print fit.a, fit.sig_a
    plot = p.errorbar([0], [fit.a], [fit.sig_a],
                      fmt='sg', lw=1, clip_on=False, alpha=0.7)
    # # Adjust error bar.
    for x in plot[1]:
        x.set_clip_on(False)
        x.set_markeredgewidth(2)

    lg = p.legend(legend, ('$am_s = 0.0102$', '$am_s=0.0509$'), 'best', numpoints=1)
    lg.get_frame().set_linewidth(0.5)
    if save:
        p.savefig(name)
    else:
        p.show()

def main(argv):
    # Load correlators.
    HHpions = [HHpion('0102', list_0102_0165), 
               HHpion('0102', list_0102_0240),
               HHpion('0102', list_0102_0380), 
               HHpion('0102', list_0102_0731),
               HHpion('0509', list_0509_0165),
               HHpion('0509', list_0509_0240),
               HHpion('0509', list_0509_0380),
               HHpion('0509', list_0509_0731)]
    
    HOpions = [HOpion('0102', '0165', list_0102_0165), 
               HOpion('0102', '0240', list_0102_0240),
               HOpion('0102', '0380', list_0102_0380),
               HOpion('0102', '0731', list_0102_0731),
               HOpion('0509', '0165', list_0509_0165),
               HOpion('0509', '0240', list_0509_0240),
               HOpion('0509', '0380', list_0509_0380),
               HOpion('0509', '0731', list_0509_0731)]

    print 'HH'
    for p in HHpions:
        print p.m1, p.m2
        f = fit_cfuns(p.JKcorrelators.real/100000, 10, 22, 64, fit_single_cosh)
        f2 = fit_cfuns_through(p.JKcorrelators.real/100000, 10, 22, 64, fit_single_cosh)
        p.fit = f
        #print p.fit
        p.msq = (p.fit[0][1])**2
        print p.msq
        p.msqJK = f2[0][:,1]**2
        print p.msqJK
        p.sig_msq = 2*p.fit[1][1] # Naive error.
    
    print 'HO'
    for p in HOpions:
        print p.m1, p.m2
        f = fit_cfuns(p.JKcorrelators.real, 15, 29, 64, fit_single_cosh_osc)
        f2 = fit_cfuns_through(p.JKcorrelators.real, 15, 29, 64, fit_single_cosh_osc)
        p.fit = f
        #print p.fit
        p.msq = (p.fit[0][2])**2
        print p.msq
        p.msqJK = f2[0][:,2]**2
        print p.msqJK
        p.sig_msq = 2*p.fit[1][2]  # Naive error.

    print 'final'
    for i in range(8):
        print HOpions[i].m1, HOpions[i].m2
        x = HOpions[i].msqJK - 0.5*HHpions[i].msqJK
        print x[0]
        print JKsigma(x)
        print ''

    plot_dmsq(HOpions, HHpions, save=True, 
                                name=sroot+'delta_msq_correlated.pdf')

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
