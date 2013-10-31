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
    
def plot_fDs(results32, results48, save=False, name='', title=None):
    x32, y32, e32 = zip(*results32) # Unpack results.
    x48, y48, e48 = zip(*results48) # Unpack results.
    
    legend = ()
    p.figure()
    if title is not None:
        p.title(title)
    p.xlabel('$mc + ms$ [MeV]')
    #p.ylabel('$f_{D_s}$ [MeV]')
    p.ylabel('$f_{D^*_s}/f_{D_s}$')
    legend += p.errorbar(x32, y32, e32, fmt='ks')[0],
    legend += p.errorbar(x48, y48,e48, fmt='bo')[0],
    p.legend(legend,('$a^{-1} = 2222$ [MeV]','$a^{-1} = 3390$ [MeV]'), 'best')
    if save:
        p.savefig(name)
    else:
        p.show()
        
def plot_fDs_single(results, save=False, name='', title=None, legend_spec=''):
    x, y, e = zip(*results) # Unpack results.
    
    legend = ()
    p.figure()
    if title is not None:
        p.title(title)
    p.xlabel('$am_s$')
    p.xlim([.026,.031])
    p.ylabel('$f_{D_s}$ [MeV]')
    legend += p.errorbar(x, y, e, fmt='ks')[0],
    #p.legend(legend, legend_spec, 'best')
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

    fDs_32_pp = []   # Aggregate results.
    ratio_32_pp = []  # Aggregate results.  
    
    for mc, ms in itertools.product(charm_masses, strange_masses):
    #for mc, ms in [(0.425, 0.048), (.430, .048), (0.435, 0.048)]:
        print '__mc = {0}, ms = {1}__'.format(mc, ms)
        x = Overlap_T96(mc, ms)
        spec = 'L32T96_c{0}_s{1}'.format(mc, ms)
        
#        print 'Point-Point pseudoscalar:'
#        name = root + spec + '_pscalar_pp_'.format(mc, ms)
#        if options.plot:
#            plot_correlator(x.pscalar_pp, options.save, name+'corr.pdf', spec)
#            plot_effmass(x.pscalar_pp, options.save, name+'meff_naive.pdf', 
#                         spec, [0.6,1.0])
#        p1, err, chisq = fit_twopoint_cfuns(x.pscalar_pp.real, 20, 42, x.T)
#        print 'A=', p1[0], '+/-', err[0]
#        print 'm=', p1[1], '+/-', err[1], '-->', p1[1]*ainv/csq, '+/-',\
#                                                 err[1]*ainv/csq
#        print 'chisq=', chisq
#        # f_Ds.
#        JK = fit_twopoint_cfunsJK(x.pscalar_pp.real, 20, 42, x.T)
#        A_JK, M_JK = np.array(zip(*JK))
#        tmpJK = f_P(mc, ms, M_JK, A_JK)
#        ave, sig = tmpJK[0], JKsigma(tmpJK)
#        print 'f_Ds = ', ave, '+/-', sig, '-->',\
#                         ave*ainv, '+/-', sig*ainv,  '\n'
#        fDs_32_pp.append((mc, ave*ainv, sig*ainv))
#       
#        print 'Point-Wall pseudoscalar:'
#        name = root + spec + '_pscalar_pw_'.format(mc, ms)
#        if options.plot:
#            plot_correlator(x.pscalar_pw, options.save, name+'corr.pdf', spec)
#            plot_effmass(x.pscalar_pw, options.save, name+'meff_naive.pdf', 
#                         spec, [0.4,1.2])
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
#            plot_effmass(x.pscalar_wp, options.save, name+'meff_naive.pdf', 
#                         spec, [0.6,1.0])
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
#            plot_effmass(x.pscalar_ww, options.save, name+'meff_naive.pdf', 
#                         spec, [0.4,1.2])
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
    
        print 'Point-Point vector:'
        name = root + spec + '_vector_pp_'.format(mc, ms)
        p1, err, chisq = fit_twopoint_cfuns(x.vector_pp.real, 15, 45, x.T)
        print 'A=', p1[0], '+/-', err[0]
        print 'm=', p1[1], '+/-', err[1], '-->', p1[1]*ainv/csq, '+/-',\
                                                 err[1]*ainv/csq
        print 'chisq=', chisq, '\n'
        if options.plot:
            plot_correlator(x.vector_pp, options.save, name+'corr.pdf', spec)
            plot_effmass(x.vector_pp, options.save, name+'meff_naive.pdf', 
                         spec, [0.6,1.2])
        
#        print 'Point-Wall vector:'
#        name = root + spec + '_vector_pw_'.format(mc, ms)
#        if options.plot:
#            plot_correlator(x.vector_pw, options.save, name+'corr.pdf', spec)
#            plot_effmass(x.vector_pw, options.save, name+'meff_naive.pdf', 
#                         spec, [0.4,1.4])
#        p1, err, chisq = fit_twopoint_cfuns(x.vector_pw.real, 15, 45, x.T)
#        print 'A=', p1[0], '+/-', err[0]
#        print 'm=', p1[1], '+/-', err[1], '-->', p1[1]*ainv/csq, '+/-',\
#                                                 err[1]*ainv/csq
#        print 'chisq=', chisq, '\n'
#             
#        print 'Wall-Point vector:'
#        name = root + spec + '_vector_wp_'.format(mc, ms)
#        if options.plot:
#            plot_correlator(x.vector_wp, options.save, name+'corr.pdf', spec)
#            plot_effmass(x.vector_wp, options.save, name+'meff_naive.pdf', 
#                         spec, [0.6,1.2])
#        p1, err, chisq = fit_twopoint_cfuns(x.vector_wp.real, 15, 45, x.T)
#        print 'A=', p1[0], '+/-', err[0]
#        print 'm=', p1[1], '+/-', err[1], '-->', p1[1]*ainv/csq, '+/-',\
#                                                 err[1]*ainv/csq
#        print 'chisq=', chisq, '\n'
#        
#        print 'Wall-Wall vector:'
#        name = root + spec + '_vector_ww_'.format(mc, ms)
#        if options.plot:
#            plot_correlator(x.vector_ww, options.save, name+'corr.pdf', spec)
#            plot_effmass(x.vector_ww, options.save, name+'meff_naive.pdf', 
#                         spec, [0.4,1.4])
#        p1, err, chisq = fit_twopoint_cfuns(x.vector_ww.real, 15, 45, x.T)
#        print 'A=', p1[0], '+/-', err[0]
#        print 'm=', p1[1], '+/-', err[1], '-->', p1[1]*ainv/csq, '+/-',\
#                                                 err[1]*ainv/csq
#        print 'chisq=', chisq, '\n'

        print 'Point-Point a4a4:'
        name = root + spec + '_a4a4_pp_'
        p1, err, chisq = fit_twopoint_cfuns(-x.a4a4_pp.real, 15, 45, x.T)
        print 'A=', p1[0], '+/-', err[0]
        print 'm=', p1[1], '+/-', err[1], '-->', p1[1]*ainv/csq, '+/-',\
                                                 err[1]*ainv/csq

        print 'chisq=', chisq, '\n'
        if options.plot:
            plot_correlator(-x.a4a4_pp, options.save, name+'corr.pdf', spec)
            plot_effmass(-x.a4a4_pp, options.save, name+'meff_naive.pdf', 
                         spec, [0.6,1.0])
                    
#        print 'Point-Wall a4a4:'
#        name = root + spec + '_a4a4_pw_'.format(mc, ms)
#        p1, err, chisq = fit_twopoint_cfuns(-x.a4a4_pw.real, 15, 45, x.T)
#        print 'A=', p1[0], '+/-', err[0]
#        print 'm=', p1[1], '+/-', err[1], '-->', p1[1]*ainv/csq, '+/-',\
#                                                 err[1]*ainv/csq
#        print 'chisq=', chisq, '\n'
#        if options.plot:
#            plot_correlator(-x.a4a4_pw, options.save, name+'corr.pdf', spec)
#            plot_effmass(-x.a4a4_pw, options.save, name+'meff_naive.pdf', 
#                         spec, [0.4,1.2])    
#         
#        print 'Wall-Point a4a4:'
#        name = root + spec + '_a4a4_wp_'.format(mc, ms)
#        p1, err, chisq = fit_twopoint_cfuns(-x.a4a4_wp.real, 15, 45, x.T)
#        print 'A=', p1[0], '+/-', err[0]
#        print 'm=', p1[1], '+/-', err[1], '-->', p1[1]*ainv/csq, '+/-',\
#                                                 err[1]*ainv/csq
#        print 'chisq=', chisq, '\n'
#        if options.plot:
#            plot_correlator(-x.a4a4_wp, options.save, name+'corr.pdf', spec)
#            plot_effmass(-x.a4a4_wp, options.save, name+'meff_naive.pdf', 
#                         spec, [0.6,1.0])
#        
#        print 'Wall-Wall a4a4:'
#        name = root + spec + '_a4a4_ww_'.format(mc, ms)
#        p1, err, chisq = fit_twopoint_cfuns(-x.a4a4_ww.real, 15, 45, x.T)
#        print 'A=', p1[0], '+/-', err[0]
#        print 'm=', p1[1], '+/-', err[1], '-->', p1[1]*ainv/csq, '+/-',\
#                                                 err[1]*ainv/csq
#        print 'chisq=', chisq, '\n'
#        if options.plot:
#            plot_correlator(-x.a4a4_ww, options.save, name+'corr.pdf', spec)
#            plot_effmass(-x.a4a4_ww, options.save, name+'meff_naive.pdf', 
#                         spec, [0.4,1.2])
        
        print 'Working on f_Ds*/f_Ds ratio...'
        r_corr = x.vector_pp.real/(-x.a4a4_pp.real)
        name = root + spec + 'vva4a4_ratio_pp'.format(mc, ms)
        #print x.vector_pp.real[0]
        #print -x.a4a4_pp.real[0]
        #print r_corr.real[0]
        p1, err, chisq = fit_twopoint_cfuns(r_corr.real, 15, 45, x.T)
        print 'A=', p1[0], '+/-', err[0]
        print 'm=', p1[1], '+/-', err[1], '-->', p1[1]*ainv/csq, '+/-',\
                                                 err[1]*ainv/csq
        print 'chisq=', chisq, '\n'
        if options.plot:
            plot_correlator(r_corr, options.save, name+'corr.pdf', spec)
            plot_effmass(r_corr, options.save, name+'meff_naive.pdf', spec)
            
        print 'Fit V-V and A4-A4 separately:'
        vector_ppJK = fit_twopoint_cfunsJK(x.vector_pp.real, 15, 45, x.T)
        vectorA_JK, vectorM_JK = np.array(zip(*vector_ppJK))
        a4a4_ppJK = fit_twopoint_cfunsJK(-x.a4a4_pp.real, 15, 45, x.T)
        a4a4A_JK, a4a4M_JK = np.array(zip(*a4a4_ppJK))
        tmpJK = f_ratio(vectorM_JK, vectorA_JK, a4a4M_JK, a4a4A_JK)
        ave, sigma = tmpJK[0], JKsigma(tmpJK)
        print 'ratio:', ave, '+/-', sigma, '\n'
        ratio_32_pp.append((ainv*(mc+ms), ave, sigma))
            
        print ''

    print np.array(ratio_32_pp)        
#        print fDs_32_pp
#    plot_fDs_single(fDs_32_pp, save=True, 
#                    name='/Users/atlytle/Desktop/fDs_s{0}.pdf'.format(ms), 
#                    title='$am_s={0}$, $1/a = 2222$ MeV'.format(ms),
#                    legend_spec=('$am_c=$'.format(mc),))
    #L48T144.
    print '------------- L48T144 -------------'
    
    ainv = 3390.48  # (23) MeV
    csq = 0.93  # c^2.
    charm_masses = [0.28, 0.29, 0.3]
    strange_masses = [0.027, 0.028, 0.029, 0.030]
    
    fDs_48_pp = []  # Aggregate results.
    fDs_wp = []
    ratio_48_pp = []  # Aggregate results.  
    for mc, ms in itertools.product(charm_masses, strange_masses):
    #for mc, ms in [(0.29, 0.027), (0.29, 0.028), (0.29, 0.029), (0.29, 0.03)]:
        print '\n__mc = {0}, ms = {1}__'.format(mc, ms)
        x = Overlap_T144(mc, ms)
        spec = 'L48T144_c{0}_s{1}'.format(mc, ms)
#        
#        print 'Point-Point pseudoscalar:'
#        name = root + spec + '_pscalar_pp_'.format(mc, ms)
#        if options.plot:
#            plot_correlator(x.pscalar_pp, options.save, name+'corr.pdf', spec)
#            plot_effmass(x.pscalar_pp, options.save, name+'meff_naive.pdf', 
#                         spec, [0.4,0.7])
#        p1, err, chisq = fit_twopoint_cfuns(x.pscalar_pp.real, 20, 42, x.T)
#        print 'A=', p1[0], '+/-', err[0]
#        print 'm=', p1[1], '+/-', err[1], '-->', p1[1]*ainv/csq, '+/-',\
#                                                 err[1]*ainv/csq
#        print 'chisq=', chisq
#        
#        # f_Ds.
#        JK = fit_twopoint_cfunsJK(x.pscalar_pp.real, 20, 42, x.T)
#        A_JK, M_JK = np.array(zip(*JK))
#        tmpJK = f_P(mc, ms, M_JK, A_JK)
#        ave, sig = tmpJK[0], JKsigma(tmpJK)
#        print 'f_Ds = ', ave, '+/-', sig, '-->',\
#                         ave*ainv, '+/-', sig*ainv,  '\n'
#        fDs_48_pp.append((ms, ave*ainv, sig*ainv))
#        
#        print 'Point-Wall pseudoscalar:'
#        name = root + spec + '_pscalar_pw_'.format(mc, ms)
#        if options.plot:
#            plot_correlator(x.pscalar_pw, options.save, name+'corr.pdf', spec)
#            plot_effmass(x.pscalar_pw, options.save, name+'meff_naive.pdf', 
#                         spec, [0.2,0.8])
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
#            plot_effmass(x.pscalar_wp, options.save, name+'meff_naive.pdf', 
#                         spec, [0.4,0.7])
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
#            plot_effmass(x.pscalar_ww, options.save, name+'meff_naive.pdf', 
#                         spec, [0.2,0.8])
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

        print 'Point-Point vector:'
        name = root + spec + '_vector_pp_'.format(mc, ms)
        if options.plot:
            plot_correlator(x.vector_pp, options.save, name+'corr.pdf', spec)
            plot_effmass(x.vector_pp, options.save, name+'meff_naive.pdf', 
                         spec, [0.2,0.8])
        p1, err, chisq = fit_twopoint_cfuns(x.vector_pp.real, 30, 55, x.T)
        print 'A=', p1[0], '+/-', err[0]
        print 'm=', p1[1], '+/-', err[1], '-->', p1[1]*ainv/csq, '+/-',\
                                                 err[1]*ainv/csq
        print 'chisq=', chisq, '\n'
#        
#        print 'Point-Wall vector:'
#        name = root + spec + '_vector_pw_'.format(mc, ms)
#        if options.plot:
#            plot_correlator(x.vector_pw, options.save, name+'corr.pdf', spec)
#            plot_effmass(x.vector_pw, options.save, name+'meff_naive.pdf', 
#                         spec, [0.2,0.8])
#        p1, err, chisq = fit_twopoint_cfuns(x.vector_pw.real, 30, 55, x.T)
#        print 'A=', p1[0], '+/-', err[0]
#        print 'm=', p1[1], '+/-', err[1], '-->', p1[1]*ainv/csq, '+/-',\
#                                                 err[1]*ainv/csq
#        print 'chisq=', chisq, '\n'
#             
#        print 'Wall-Point vector:'
#        name = root + spec + '_vector_wp_'.format(mc, ms)
#        if options.plot:
#            plot_correlator(x.vector_wp, options.save, name+'corr.pdf', spec)
#            plot_effmass(x.vector_wp, options.save, name+'meff_naive.pdf', 
#                         spec, [0.2,0.8])
#        p1, err, chisq = fit_twopoint_cfuns(x.vector_wp.real, 30, 55, x.T)
#        print 'A=', p1[0], '+/-', err[0]
#        print 'm=', p1[1], '+/-', err[1], '-->', p1[1]*ainv/csq, '+/-',\
#                                                 err[1]*ainv/csq
#        print 'chisq=', chisq, '\n'
#        
#        print 'Wall-Wall vector:'
#        name = root + spec + '_vector_ww_'.format(mc, ms)
#        if options.plot:
#            plot_correlator(x.vector_ww, options.save, name+'corr.pdf', spec)
#            plot_effmass(x.vector_ww, options.save, name+'meff_naive.pdf', 
#                         spec, [0.2,0.8])
#        p1, err, chisq = fit_twopoint_cfuns(x.vector_ww.real, 30, 55, x.T)
#        print 'A=', p1[0], '+/-', err[0]
#        print 'm=', p1[1], '+/-', err[1], '-->', p1[1]*ainv/csq, '+/-',\
#                                                 err[1]*ainv/csq
#        print 'chisq=', chisq, '\n'


        print 'Point-Point a4a4:'
        name = root + spec + '_a4a4_pp_'
        if options.plot:
            plot_correlator(-x.a4a4_pp, options.save, name+'corr.pdf', spec)
            plot_effmass(-x.a4a4_pp, options.save, name+'meff_naive.pdf', 
                         spec, [0.4,0.7])
        p1, err, chisq = fit_twopoint_cfuns(-x.a4a4_pp.real, 30, 55, x.T)
        print 'A=', p1[0], '+/-', err[0]
        print 'm=', p1[1], '+/-', err[1], '-->', p1[1]*ainv/csq, '+/-',\
                                                 err[1]*ainv/csq
        print 'chisq=', chisq, '\n'
#        
#        print 'Point-Wall a4a4:'
#        name = root + spec + '_a4a4_pw_'.format(mc, ms)
#        if options.plot:
#            plot_correlator(-x.a4a4_pw, options.save, name+'corr.pdf', spec)
#            plot_effmass(-x.a4a4_pw, options.save, name+'meff_naive.pdf', 
#                         spec, [0.2,0.8])
#        p1, err, chisq = fit_twopoint_cfuns(-x.a4a4_pw.real, 30, 55, x.T)
#        print 'A=', p1[0], '+/-', err[0]
#        print 'm=', p1[1], '+/-', err[1], '-->', p1[1]*ainv/csq, '+/-',\
#                                                 err[1]*ainv/csq
#        print 'chisq=', chisq, '\n'
#             
#        print 'Wall-Point a4a4:'
#        name = root + spec + '_a4a4_wp_'.format(mc, ms)
#        if options.plot:
#            plot_correlator(-x.a4a4_wp, options.save, name+'corr.pdf', spec)
#            plot_effmass(-x.a4a4_wp, options.save, name+'meff_naive.pdf', 
#                         spec, [0.4,0.7])
#        p1, err, chisq = fit_twopoint_cfuns(-x.a4a4_wp.real, 30, 55, x.T)
#        print 'A=', p1[0], '+/-', err[0]
#        print 'm=', p1[1], '+/-', err[1], '-->', p1[1]*ainv/csq, '+/-',\
#                                                 err[1]*ainv/csq
#        print 'chisq=', chisq, '\n'
#        
#        print 'Wall-Wall a4a4:'
#        name = root + spec + '_a4a4_ww_'.format(mc, ms)
#        if options.plot:
#            plot_correlator(-x.a4a4_ww, options.save, name+'corr.pdf', spec)
#            plot_effmass(-x.a4a4_ww, options.save, name+'meff_naive.pdf', 
#                         spec, [0.2,0.8])
#        p1, err, chisq = fit_twopoint_cfuns(-x.a4a4_ww.real, 30, 55, x.T)
#        print 'A=', p1[0], '+/-', err[0]
#        print 'm=', p1[1], '+/-', err[1], '-->', p1[1]*ainv/csq, '+/-',\
#                                                 err[1]*ainv/csq
#        print 'chisq=', chisq, '\n'

        print 'Working on f_Ds*/f_Ds ratio...'
        r_corr = x.vector_pp.real/(-x.a4a4_pp.real)
        name = root + spec + 'vva4a4_ratio_pp'.format(mc, ms)
        #print x.vector_pp.real[0]
        #print -x.a4a4_pp.real[0]
        #print r_corr.real[0]
        if options.plot:
            plot_correlator(r_corr, options.save, name+'corr.pdf', spec)
            plot_effmass(r_corr, options.save, name+'meff_naive.pdf', spec)
        p1, err, chisq = fit_twopoint_cfuns(r_corr.real, 30, 55, x.T)
        print 'A=', p1[0], '+/-', err[0]
        print 'm=', p1[1], '+/-', err[1], '-->', p1[1]*ainv/csq, '+/-',\
                                                 err[1]*ainv/csq
        print 'chisq=', chisq, '\n'
        
        print 'Fit V-V and A4-A4 separately:'
        vector_ppJK = fit_twopoint_cfunsJK(x.vector_pp.real, 30, 55, x.T)
        vectorA_JK, vectorM_JK = np.array(zip(*vector_ppJK))
        a4a4_ppJK = fit_twopoint_cfunsJK(-x.a4a4_pp.real, 30, 55, x.T)
        a4a4A_JK, a4a4M_JK = np.array(zip(*a4a4_ppJK))
        tmpJK = f_ratio(vectorM_JK, vectorA_JK, a4a4M_JK, a4a4A_JK)
        ave, sigma = tmpJK[0], JKsigma(tmpJK)
        print 'ratio:', ave, '+/-', sigma, '\n'
        ratio_48_pp.append((ainv*(mc+ms), ave, sigma))
        
        print ''
    print np.array(ratio_48_pp)
    plot_fDs(ratio_32_pp, ratio_48_pp, True, root+'fratio_pp.pdf')
#    print fDs_48_pp
#    plot_fDs(fDs_32_pp, fDs_48_pp, True, root+'fDs_pscalar_pp.pdf')
#    plot_fDs_single(fDs_48_pp, save=False, 
#                    name='/Users/atlytle/Desktop/fDs_c{0}.pdf'.format(ms), 
#                    title='$am_c={0}$, $1/a = 3390$ MeV'.format(mc)) 
    return 0
    
if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
