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
    
def plot_fDs(results32, results48, save=False, name='', title=None):
    x32, y32, e32 = zip(*results32) # Unpack results.
    x48, y48, e48 = zip(*results48) # Unpack results.
    
    legend = ()
    p.figure()
    if title is not None:
        p.title(title)
    p.xlabel('$mc + ms$ [MeV]')
    p.ylabel('$f_{D_s}$ [MeV]')
    legend += p.errorbar(x32, y32, e32, fmt='ks')[0],
    legend += p.errorbar(x48, y48,e48, fmt='bo')[0],
    p.legend(legend,('$a^{-1} = 2222$ [MeV]','$a^{-1} = 3390$ [MeV]'), 'best')
    if save:
        p.savefig(name)
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

    fDs_32_pp = []   #Aggregate results.
    
    for mc, ms in itertools.product(charm_masses, strange_masses):
        print '__mc = {0}, ms = {1}__'.format(mc, ms)
        x = Overlap_T96(mc, ms)
        print x.N
        spec = 'L32T96_c{0}_s{1}'.format(mc, ms)
        
        print 'Point-Point pseudoscalar:'
        name = root + spec + '_pscalar_pp_'.format(mc, ms)
        if options.plot:
            plot_correlator(x.pscalar_pp, options.save, name+'corr.pdf', spec)
            plot_effmass(x.pscalar_pp, options.save, name+'meff_naive.pdf', spec)
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
        fDs_32_pp.append((ainv*(mc+ms), ave*ainv, sig*ainv))
       
#        print 'Point-Wall pseudoscalar:'
#        name = root + spec + '_pscalar_pw_'.format(mc, ms)
#        if options.plot:
#            plot_correlator(x.pscalar_pw, options.save, name+'corr.pdf', spec)
#            plot_effmass(x.pscalar_pw, options.save, name+'meff_naive.pdf', spec)
#        p1, err, chisq = fit_twopoint_cfuns(x.pscalar_pw.real, 20, 42, x.T)
#        print 'A=', p1[0], '+/-', err[0]
#        print 'm=', p1[1], '+/-', err[1], '-->', p1[1]*ainv/csq, '+/-',\
#                                                 err[1]*ainv/csq
#        print 'chisq=', chisq, '\n'
#        
#        print 'Wall-Point pseudoscalar:'
#        name = root + spec + '_pscalar_wp_'.format(mc, ms)
#        if options.plot:
#            plot_correlator(x.pscalar_wp, options.save, name+'corr.pdf', spec)
#            plot_effmass(x.pscalar_wp, options.save, name+'meff_naive.pdf', spec)
#        p1, err, chisq = fit_twopoint_cfuns(x.pscalar_wp.real, 20, 42, x.T)
#        x.pscalar_wp_A = p1[0], err[0]
#        x.pscalar_wp_m = p1[1], err[1]
#        print 'A=', p1[0], '+/-', err[0]
#        print 'm=', p1[1], '+/-', err[1], '-->', p1[1]*ainv/csq, '+/-',\
#                                                 err[1]*ainv/csq
#        print 'chisq=', chisq, '\n'
#        
#        print 'Wall-Wall pseudoscalar:'
#        name = root + spec + '_pscalar_ww_'.format(mc, ms)
#        if options.plot:
#            plot_correlator(x.pscalar_ww, options.save, name+'corr.pdf', spec)
#            plot_effmass(x.pscalar_ww, options.save, name+'meff_naive.pdf', spec)
#        p1, err, chisq = fit_twopoint_cfuns(x.pscalar_ww.real, 20, 42, x.T)
#        x.pscalar_ww_A = p1[0], err[0]
#        x.pscalar_ww_m = p1[1], err[1]
#        print 'A=', p1[0], '+/-', err[0]
#        print 'm=', p1[1], '+/-', err[1], '-->', p1[1]*ainv/csq, '+/-',\
#                                                 err[1]*ainv/csq
#        print 'chisq=', chisq
#        
#        Apw, Aww = x.pscalar_wp_A[0], x.pscalar_ww_A[0]
#        tmp = Apw*Apw/(Aww*32*32*32)
#        print 'A(pp) via A(wp) and A(ww):', tmp 
#        tmp2 = f_P(mc, ms, x.pscalar_wp_m[0],tmp)
#        print 'f_Ds via A(wp) and A(ww):', tmp2, '-->', tmp2*ainv
#        print ''
#    print fDs_32_pp
#    #plot_fDs(fDs_32_pp)
        
    # L48T144.
    print '------------- L48T144 -------------'
    
    ainv = 3390.48  # (23) MeV
    csq = 0.93  # c^2.
    charm_masses = [0.28, 0.29, 0.3]
    strange_masses = [0.027, 0.028, 0.029, 0.030]
    
    fDs_48_pp = []  # Aggregate results.
    fDs_wp = []
    for mc, ms in itertools.product(charm_masses, strange_masses):
        print '\n__mc = {0}, ms = {1}__'.format(mc, ms)
        x = Overlap_T144(mc, ms)
        spec = 'L48T144_c{0}_s{1}'.format(mc, ms)
        
        print 'Point-Point pseudoscalar:'
        name = root + spec + '_pscalar_pp_'.format(mc, ms)
        if options.plot:
            plot_correlator(x.pscalar_pp, options.save, name+'corr.pdf', spec)
            plot_effmass(x.pscalar_pp, options.save, name+'meff_naive.pdf', spec)
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
        fDs_48_pp.append((ainv*(mc+ms), ave*ainv, sig*ainv))
#        
#        print 'Point-Wall pseudoscalar:'
#        name = root + spec + '_pscalar_pw_'.format(mc, ms)
#        if options.plot:
#            plot_correlator(x.pscalar_pw, options.save, name+'corr.pdf', spec)
#            plot_effmass(x.pscalar_pw, options.save, name+'meff_naive.pdf', spec)
#        p1, err, chisq = fit_twopoint_cfuns(x.pscalar_pw.real, 20, 42, x.T)
#        print 'A=', p1[0], '+/-', err[0]
#        print 'm=', p1[1], '+/-', err[1], '-->', p1[1]*ainv/csq, '+/-',\
#                                                 err[1]*ainv/csq
#        print 'chisq=', chisq, '\n'
#        
#        print 'Wall-Point pseudoscalar:'
#        name = root + spec + '_pscalar_wp_'.format(mc, ms)
#        if options.plot:
#            plot_correlator(x.pscalar_wp, options.save, name+'corr.pdf', spec)
#            plot_effmass(x.pscalar_wp, options.save, name+'meff_naive.pdf', spec)
#        p1, err, chisq = fit_twopoint_cfuns(x.pscalar_wp.real, 20, 42, x.T)
#        x.pscalar_wp_A = p1[0], err[0]
#        x.pscalar_wp_m = p1[1], err[1]
#        print 'A=', p1[0], '+/-', err[0]
#        print 'm=', p1[1], '+/-', err[1], '-->', p1[1]*ainv/csq, '+/-',\
#                                                 err[1]*ainv/csq
#        print 'chisq=', chisq, '\n'
#        
#        print 'Wall-Wall pseudoscalar:'
#        name = root + spec + '_pscalar_ww_'.format(mc, ms)
#        if options.plot:
#            plot_correlator(x.pscalar_ww, options.save, name+'corr.pdf', spec)
#            plot_effmass(x.pscalar_ww, options.save, name+'meff_naive.pdf', spec)
#        p1, err, chisq = fit_twopoint_cfuns(x.pscalar_ww.real, 30, 70, x.T)
#        x.pscalar_ww_A = p1[0], err[0]
#        x.pscalar_ww_m = p1[1], err[1]
#        print 'A=', p1[0], '+/-', err[0]
#        print 'm=', p1[1], '+/-', err[1], '-->', p1[1]*ainv/csq, '+/-',\
#                                                 err[1]*ainv/csq
#        print 'chisq=', chisq
#        
#        Apw, Aww = x.pscalar_wp_A[0], x.pscalar_ww_A[0]
#        tmp = Apw*Apw/(Aww*48*48*48)
#        print 'A(pp) via A(wp) and A(ww):', tmp 
#        tmp2 = f_P(mc, ms, x.pscalar_wp_m[0],tmp)
#        print 'f_Ds via A(wp) and A(ww):', tmp2, '-->', tmp2*ainv
#        fDs_wp.append((mc+ms, tmp2*ainv))
#        print ''
#    print fDs_48_pp
    plot_fDs(fDs_32_pp, fDs_48_pp, True, root+'fDs_pscalar_pp.pdf')
        
    return 0
    
if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
