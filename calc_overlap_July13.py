import optparse
import sys
import itertools

import numpy as np

from parse_overlap import Overlap_T96, Overlap_T144
from calc_overlap import plot_correlator, plot_effmass, fit_twopoint_cfuns

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
    
def main(argv=None):
    options = parse_args(argv)
    root = '/Users/atlytle/Dropbox/TIFR/figs/'
    
#    # L32T96.
#    print '---------- L32T96 ------------'
#    ainv = 2222.15  # (20) MeV.
#    csq = 0.85  # c^2.
#    charm_masses = [0.425, 0.430, 0.435]
#    strange_masses = [0.046, 0.048, 0.049]
#    
#    for mc, ms in itertools.product(charm_masses, strange_masses):
#        print 'mc = {0}, ms = {1}:'.format(mc, ms)
#        x = Overlap_T96(mc, ms)
#        spec = 'L32T96_c{0}_s{1}'.format(mc, ms)
#        
#        print 'Point-Point pseudoscalar:'
#        name = root + spec + '_pscalar_pp_'.format(mc, ms)
#        if options.plot:
#            plot_correlator(x.pscalar_pp, options.save, name+'corr.pdf', spec)
#            plot_effmass(x.pscalar_pp, options.save, name+'meff_naive.pdf', spec)
#        p1, err, chisq = fit_twopoint_cfuns(x.pscalar_pp.real, 20, 42, x.T)
#        print 'A=', p1[0], '+/-', err[0]
#        print 'm=', p1[1], '+/-', err[1], '-->', p1[1]*ainv/csq, '+/-',\
#                                                 err[1]*ainv/csq
#        print 'chisq=', chisq
#        tmp = f_P(mc, ms, p1[1], p1[0])
#        print 'f_Ds = ', tmp, '-->', tmp*ainv,  '\n'
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

        
    # L48T144.
    print '------------- L48T144 -------------'
    
    ainv = 3390.48  # (23) MeV
    csq = 0.93  # c^2.
    charm_masses = [0.28, 0.29, 0.3]
    strange_masses = [0.027, 0.028, 0.029, 0.030]
    
    for mc, ms in itertools.product(charm_masses, strange_masses):
        print 'mc = {0}, ms = {1}:'.format(mc, ms)
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
        tmp = f_P(mc, ms, p1[1], p1[0])
        print 'f_Ds = ', tmp, '-->', tmp*ainv,  '\n'
        print 'Point-Wall pseudoscalar:'
        name = root + spec + '_pscalar_pw_'.format(mc, ms)
        if options.plot:
            plot_correlator(x.pscalar_pw, options.save, name+'corr.pdf', spec)
            plot_effmass(x.pscalar_pw, options.save, name+'meff_naive.pdf', spec)
        p1, err, chisq = fit_twopoint_cfuns(x.pscalar_pw.real, 20, 42, x.T)
        print 'A=', p1[0], '+/-', err[0]
        print 'm=', p1[1], '+/-', err[1], '-->', p1[1]*ainv/csq, '+/-',\
                                                 err[1]*ainv/csq
        print 'chisq=', chisq, '\n'
        
        print 'Wall-Point pseudoscalar:'
        name = root + spec + '_pscalar_wp_'.format(mc, ms)
        if options.plot:
            plot_correlator(x.pscalar_wp, options.save, name+'corr.pdf', spec)
            plot_effmass(x.pscalar_wp, options.save, name+'meff_naive.pdf', spec)
        p1, err, chisq = fit_twopoint_cfuns(x.pscalar_wp.real, 20, 42, x.T)
        x.pscalar_wp_A = p1[0], err[0]
        x.pscalar_wp_m = p1[1], err[1]
        print 'A=', p1[0], '+/-', err[0]
        print 'm=', p1[1], '+/-', err[1], '-->', p1[1]*ainv/csq, '+/-',\
                                                 err[1]*ainv/csq
        print 'chisq=', chisq, '\n'
        
        print 'Wall-Wall pseudoscalar:'
        name = root + spec + '_pscalar_ww_'.format(mc, ms)
        if options.plot:
            plot_correlator(x.pscalar_ww, options.save, name+'corr.pdf', spec)
            plot_effmass(x.pscalar_ww, options.save, name+'meff_naive.pdf', spec)
        p1, err, chisq = fit_twopoint_cfuns(x.pscalar_ww.real, 20, 42, x.T)
        x.pscalar_ww_A = p1[0], err[0]
        x.pscalar_ww_m = p1[1], err[1]
        print 'A=', p1[0], '+/-', err[0]
        print 'm=', p1[1], '+/-', err[1], '-->', p1[1]*ainv/csq, '+/-',\
                                                 err[1]*ainv/csq
        print 'chisq=', chisq
        
        Apw, Aww = x.pscalar_wp_A[0], x.pscalar_ww_A[0]
        tmp = Apw*Apw/(Aww*48*48*48)
        print 'A(pp) via A(wp) and A(ww):', tmp 
        tmp2 = f_P(mc, ms, x.pscalar_wp_m[0],tmp)
        print 'f_Ds via A(wp) and A(ww):', tmp2, '-->', tmp2*ainv
        print ''
        
    return 0
    
if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
