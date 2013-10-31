import optparse
import sys
import itertools
import pylab as p

import numpy as np

from parse_overlap import Overlap_T96, Overlap_T144
from calc_overlap import plot_correlator, plot_effmass, fit_twopoint_cfuns,\
                         fit_twopoint_cfunsJK, JKsigma

def f_P(mc, ms, mP, A):
    'Pseudoscalar decay constant from |<0|P|P>| matrix element.'
    return ((mc+ms)/mP**2)*np.sqrt(2*A*mP)
    
def f_ratio(mV, AV, mP, AP):
    'fDs*/fDs ratio from V-V and A4-A4 correlators. Ignores ZA/ZV.'
    tmp = ((AV/3)/AP)*(mP/mV)
    return np.sqrt(tmp)
    
def parse_args(argv):
    parser = optparse.OptionParser()
    parser.add_option('-p', '--plot', action='store_true', dest='plot',
                      help = 'Plot results.')
    parser.add_option('-s', '--save', action='store_true', dest='save',
                      help = 'Save the plots.')
    options, args = parser.parse_args(argv)
    return options
    
def filestring(mc, ms, spec):
    return 'c{0}_u{1}_{2}.pdf'.format(mc, ms, spec)
        
def plot_fDs(results, save=False, saveas='', title=None, legend_spec='',
             xlabel='', xlim=[], ylim=[]):
    
    x, y, e = zip(*results) # Unpack results.
    
    legend = ()
    p.figure()
    p.rc('text', usetex=True)
    p.rc('font', size=18)
    if title is not None:
        p.title(title)
    if xlabel:
        p.xlabel(xlabel)
    if xlim:
        p.xlim(xlim)
    p.ylabel('$f_{D_s} \mathrm{\,\, [MeV]}$')
    if ylim:
        p.ylim(ylim)
    p.errorbar(x, y, e, fmt='rs', ms=10),
    p.axhspan((248.6-2.7), (248.6+2.7),facecolor=(0.,0.5,1.), alpha=.3)
    if save:
        p.savefig(saveas)
    else:
        p.show()
        
def plot_fratio(results, save=False, saveas='', title=None, legend_spec='',
             xlabel='', xlim=[], ylim=[]):
    
    x, y, e = zip(*results) # Unpack results.
    
    legend = ()
    p.figure()
    if title is not None:
        p.title(title)
    if xlabel:
        p.xlabel(xlabel)
    if xlim:
        p.xlim(xlim)
    p.ylabel('$f_{D^*_s}/f_{D_s}$')
    if ylim:
        p.ylim(ylim)
    p.errorbar(x, y, e, fmt='bo', ms=10),
    if save:
        p.savefig(saveas)
    else:
        p.show()
    
def main(argv=None):
    options = parse_args(argv)
    root = '/Users/atlytle/Dropbox/TIFR/figs/'
        
    # L32T96.
    print '---------- L32T96 ------------'
    ainv = 2222.15  # (20) MeV.
    csq = 0.85  # c^2.
    charm_masses = [0.425, 0.430, 0.435]
    strange_masses = [0.046, 0.048, 0.049]
    
    # mc = 0.430
    fDs = []  # Aggregate results.
    ratio = []
    for mc, ms in [(0.430, 0.046), (0.430, .048), (0.430, 0.049)]:
        print '__mc = {0}, ms = {1}__'.format(mc, ms)
        x = Overlap_T96(mc, ms)
        spec = 'L32T96_c{0}_s{1}'.format(mc, ms)
        
        print 'Point-Point pseudoscalar:'
        name = root + spec + '_pscalar_pp_'.format(mc, ms)
        if options.plot:
            plot_correlator(x.pscalar_pp, options.save, name+'corr.pdf', spec)
            plot_effmass(x.pscalar_pp, options.save, name+'meff_naive.pdf', 
                         spec, [0.6,1.0])
        p1, err, chisq = fit_twopoint_cfuns(x.pscalar_pp.real, 20, 42, x.T)
        print 'A=', p1[0], '+/-', err[0]
        print 'm=', p1[1], '+/-', err[1], '-->', p1[1]*ainv/csq, '+/-',\
                                                 err[1]*ainv/csq
        print 'chisq=', chisq
        # f_Ds.
        JK = fit_twopoint_cfunsJK(x.pscalar_pp.real, 20, 42, x.T)
        A_JK, M_JK = np.array(zip(*JK))
        tmpJK = f_P(mc, ms, M_JK, A_JK)
        ave, sig = tmpJK[0], JKsigma(tmpJK)
        print 'f_Ds = ', ave, '+/-', sig, '-->',\
                         ave*ainv, '+/-', sig*ainv,  '\n'
        fDs.append((ms, ave*ainv, sig*ainv))
        
        # f_Ds*/f_Ds ratio.
        print 'Fit V-V and A4-A4 separately:'
        vector_ppJK = fit_twopoint_cfunsJK(x.vector_pp.real, 15, 45, x.T)
        vectorA_JK, vectorM_JK = np.array(zip(*vector_ppJK))
        a4a4_ppJK = fit_twopoint_cfunsJK(-x.a4a4_pp.real, 15, 45, x.T)
        a4a4A_JK, a4a4M_JK = np.array(zip(*a4a4_ppJK))
        tmpJK = f_ratio(vectorM_JK, vectorA_JK, a4a4M_JK, a4a4A_JK)
        ave, sigma = tmpJK[0], JKsigma(tmpJK)
        print 'ratio:', ave, '+/-', sigma, '\n'
        ratio.append((ms, ave, sigma))
        
        
        
    plot_fDs(fDs, save=True, 
             saveas='/Users/atlytle/Desktop/fDs_c{0}.pdf'.format(mc), 
             title='$am_c={0}$, $1/a = 2222 \,\, \mathrm{{MeV}}$'.format(mc), xlabel='$am_s$',
             xlim=[.045,.050], ylim=[])
             
    plot_fratio(ratio, save=True, 
                saveas='/Users/atlytle/Desktop/fratio_c{0}.pdf'.format(mc), 
                title='$am_c={0}$, $1/a = 2222 \,\, \mathrm{{MeV}}$'.format(mc), xlabel='$am_s$',
                xlim=[.045,.050], ylim=[0.9,1.5])

    # ms = 0.048
    fDs = []  # Aggregate results.
    ratio = [] # Aggregate results.
    for mc, ms in [(0.425, 0.048), (.430, .048), (0.435, 0.048)]:
        print '__mc = {0}, ms = {1}__'.format(mc, ms)
        x = Overlap_T96(mc, ms)
        spec = 'L32T96_c{0}_s{1}'.format(mc, ms)
        
        print 'Point-Point pseudoscalar:'
        name = root + spec + '_pscalar_pp_'.format(mc, ms)
        if options.plot:
            plot_correlator(x.pscalar_pp, options.save, name+'corr.pdf', spec)
            plot_effmass(x.pscalar_pp, options.save, name+'meff_naive.pdf', 
                         spec, [0.6,1.0])
        p1, err, chisq = fit_twopoint_cfuns(x.pscalar_pp.real, 20, 42, x.T)
        print 'A=', p1[0], '+/-', err[0]
        print 'm=', p1[1], '+/-', err[1], '-->', p1[1]*ainv/csq, '+/-',\
                                                 err[1]*ainv/csq
        print 'chisq=', chisq
        # f_Ds.
        JK = fit_twopoint_cfunsJK(x.pscalar_pp.real, 20, 42, x.T)
        A_JK, M_JK = np.array(zip(*JK))
        tmpJK = f_P(mc, ms, M_JK, A_JK)
        ave, sig = tmpJK[0], JKsigma(tmpJK)
        print 'f_Ds = ', ave, '+/-', sig, '-->',\
                         ave*ainv, '+/-', sig*ainv,  '\n'
        fDs.append((mc, ave*ainv, sig*ainv))
        
        # f_Ds*/f_Ds ratio.
        print 'Fit V-V and A4-A4 separately:'
        vector_ppJK = fit_twopoint_cfunsJK(x.vector_pp.real, 15, 45, x.T)
        vectorA_JK, vectorM_JK = np.array(zip(*vector_ppJK))
        a4a4_ppJK = fit_twopoint_cfunsJK(-x.a4a4_pp.real, 15, 45, x.T)
        a4a4A_JK, a4a4M_JK = np.array(zip(*a4a4_ppJK))
        tmpJK = f_ratio(vectorM_JK, vectorA_JK, a4a4M_JK, a4a4A_JK)
        ave, sigma = tmpJK[0], JKsigma(tmpJK)
        print 'ratio:', ave, '+/-', sigma, '\n'
        ratio.append((mc, ave, sigma))
        
    plot_fDs(fDs, save=True, 
             saveas='/Users/atlytle/Desktop/fDs_s{0}.pdf'.format(ms), 
             title='$am_s={0}$, $1/a = 2222 \,\, \mathrm{{MeV}}$'.format(ms), xlabel='$am_c$',
             xlim=[], ylim=[225,260])
             
    plot_fratio(ratio, save=True, 
                saveas='/Users/atlytle/Desktop/fratio_s{0}.pdf'.format(ms), 
                title='$am_s={0}$, $1/a = 2222 \,\, \mathrm{{MeV}}$'.format(ms), xlabel='$am_c$',
                xlim=[], ylim=[0.9,1.5])
             

    
    # L48T144.
    print '------------- L48T144 -------------'
    
    ainv = 3390.48  # (23) MeV
    csq = 0.93  # c^2.
    charm_masses = [0.28, 0.29, 0.3]
    strange_masses = [0.027, 0.028, 0.029, 0.030]
    
    # mc = 0.29
    fDs = []  # Aggregate results.
    ratio = []  # Aggregate results.
    for mc, ms in [(0.29,0.027), (0.29, 0.028), (0.29, 0.029), (0.29, 0.03)]:
        print '\n__mc = {0}, ms = {1}__'.format(mc, ms)
        x = Overlap_T144(mc, ms)
        spec = 'L48T144_c{0}_s{1}'.format(mc, ms)
        
        print 'Point-Point pseudoscalar:'
        name = root + spec + '_pscalar_pp_'.format(mc, ms)
        if options.plot:
            plot_correlator(x.pscalar_pp, options.save, name+'corr.pdf', spec)
            plot_effmass(x.pscalar_pp, options.save, name+'meff_naive.pdf', 
                         spec, [0.4,0.7])
        p1, err, chisq = fit_twopoint_cfuns(x.pscalar_pp.real, 20, 42, x.T)
        print 'A=', p1[0], '+/-', err[0]
        print 'm=', p1[1], '+/-', err[1], '-->', p1[1]*ainv/csq, '+/-',\
                                                 err[1]*ainv/csq
        print 'chisq=', chisq
        
        # f_Ds.
        JK = fit_twopoint_cfunsJK(x.pscalar_pp.real, 20, 42, x.T)
        A_JK, M_JK = np.array(zip(*JK))
        tmpJK = f_P(mc, ms, M_JK, A_JK)
        ave, sig = tmpJK[0], JKsigma(tmpJK)
        print 'f_Ds = ', ave, '+/-', sig, '-->',\
                         ave*ainv, '+/-', sig*ainv,  '\n'
        fDs.append((ms, ave*ainv, sig*ainv))
        
        # f_Ds*/f_Ds ratio.
        print 'Fit V-V and A4-A4 separately:'
        vector_ppJK = fit_twopoint_cfunsJK(x.vector_pp.real, 30, 55, x.T)
        vectorA_JK, vectorM_JK = np.array(zip(*vector_ppJK))
        a4a4_ppJK = fit_twopoint_cfunsJK(-x.a4a4_pp.real, 30, 55, x.T)
        a4a4A_JK, a4a4M_JK = np.array(zip(*a4a4_ppJK))
        tmpJK = f_ratio(vectorM_JK, vectorA_JK, a4a4M_JK, a4a4A_JK)
        ave, sigma = tmpJK[0], JKsigma(tmpJK)
        print 'ratio:', ave, '+/-', sigma, '\n'
        ratio.append((ms, ave, sigma))
        
    plot_fDs(fDs, save=True, 
             saveas='/Users/atlytle/Desktop/fDs_c{0}.pdf'.format(mc), 
             title='$am_c={0}$, $1/a = 3390 \,\, \mathrm{{MeV}}$'.format(mc), xlabel='$am_s$',
             xlim=[.026,.031], ylim=[210,270]) 
             
    plot_fratio(ratio[1:], save=True, 
                saveas='/Users/atlytle/Desktop/fratio_c{0}.pdf'.format(mc), 
                title='$am_c={0}$, $1/a = 3390 \,\, \mathrm{{MeV}}$'.format(mc), xlabel='$am_s$',
                xlim=[.0275,.0305], ylim=[0.9,1.5]) 
             
    # ms = 0.028         
    fDs = []  # Aggregate results. 
    ratio = []      
    for mc, ms in [(0.28, 0.028), (0.29, 0.028), (0.3, 0.028)]:
        print '\n__mc = {0}, ms = {1}__'.format(mc, ms)
        x = Overlap_T144(mc, ms)
        spec = 'L48T144_c{0}_s{1}'.format(mc, ms)
        
        print 'Point-Point pseudoscalar:'
        name = root + spec + '_pscalar_pp_'.format(mc, ms)
        if options.plot:
            plot_correlator(x.pscalar_pp, options.save, name+'corr.pdf', spec)
            plot_effmass(x.pscalar_pp, options.save, name+'meff_naive.pdf', 
                         spec, [0.4,0.7])
        p1, err, chisq = fit_twopoint_cfuns(x.pscalar_pp.real, 20, 42, x.T)
        print 'A=', p1[0], '+/-', err[0]
        print 'm=', p1[1], '+/-', err[1], '-->', p1[1]*ainv/csq, '+/-',\
                                                 err[1]*ainv/csq
        print 'chisq=', chisq
        
        # f_Ds.
        JK = fit_twopoint_cfunsJK(x.pscalar_pp.real, 20, 42, x.T)
        A_JK, M_JK = np.array(zip(*JK))
        tmpJK = f_P(mc, ms, M_JK, A_JK)
        ave, sig = tmpJK[0], JKsigma(tmpJK)
        print 'f_Ds = ', ave, '+/-', sig, '-->',\
                         ave*ainv, '+/-', sig*ainv,  '\n'
        fDs.append((mc, ave*ainv, sig*ainv))
        
        # f_Ds*/f_Ds ratio.
        print 'Fit V-V and A4-A4 separately:'
        vector_ppJK = fit_twopoint_cfunsJK(x.vector_pp.real, 30, 55, x.T)
        vectorA_JK, vectorM_JK = np.array(zip(*vector_ppJK))
        a4a4_ppJK = fit_twopoint_cfunsJK(-x.a4a4_pp.real, 30, 55, x.T)
        a4a4A_JK, a4a4M_JK = np.array(zip(*a4a4_ppJK))
        tmpJK = f_ratio(vectorM_JK, vectorA_JK, a4a4M_JK, a4a4A_JK)
        ave, sigma = tmpJK[0], JKsigma(tmpJK)
        print 'ratio:', ave, '+/-', sigma, '\n'
        ratio.append((mc, ave, sigma))
        
    plot_fDs(fDs, save=True, 
             saveas='/Users/atlytle/Desktop/fDs_s{0}.pdf'.format(ms), 
             title='$am_s={0}$, $1/a = 3390 \,\, \mathrm{{MeV}}$'.format(ms), xlabel='$am_c$',
             xlim=[.275,.305], ylim=[225,260])
             
    plot_fratio(ratio, save=True, 
                saveas='/Users/atlytle/Desktop/fratio_s{0}.pdf'.format(ms), 
                title='$am_s={0}$, $1/a = 3390 \,\, \mathrm{{MeV}}$'.format(ms), xlabel='$am_c$',
                xlim=[.275,.305], ylim=[0.9,1.5])
    
    return 0
    
if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
