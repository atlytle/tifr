import sys
import optparse
import pylab as p
import numpy as np
import StaggeredExceptional as se
import staggered_NPR as sn

ar = np.array
np.set_printoptions(precision=2, suppress=True)

plist_stout = [(1,1,1,2), (2,2,2,2), (2,2,2,3), (2,2,2,4), (3,2,2,3),
               (2,3,3,3), (3,2,3,4), (3,3,3,3)]
gflist_stout = [100, 200, 300, 400]

plist_HYP = [(1,1,1,4), (1,1,1,6), (1,2,1,5), (1,2,2,4), (2,1,2,6), 
             (2,2,2,5), (2,2,2,7), (2,2,2,8), (2,2,2,9), (2,3,2,7)]
gflist_HYP = [600, 630, 660, 690]

plist_naive = plist_HYP
gflist_naive = [300, 330, 360, 390]

def load_data_stout(ensemble, mass, plist, gflist):
    return [se.StoutExceptional(ensemble, mass, p, gflist)
            for p in plist]
            
def load_data_HYP(ensemble, mass, plist, gflist):
    return [se.HYPExceptional(ensemble, mass, p, gflist)
            for p in plist]
            
def load_data_naive(ensemble, mass, plist, gflist):
    # HYP class has same L, T as naive and this is the only changing feature.
    return [se.HYPExceptional(ensemble, mass, p, gflist)
            for p in plist]
            
def parse_args():
    parser = optparse.OptionParser()
#    parser.add_option('-c', '--compute', action='store_true', dest='compute',
#                      help = 'Compute Zs from scratch.')
#    parser.add_option('-d', '--dump', action='store_true', dest='dump',
#                      help = 'Pickle results of the computation.')
#    parser.add_option('-l', '--load', action='store_true', dest='load',
#                      help = 'Load Zs from pickle directory.')
#    parser.add_option('-p', '--plot', action='store_true', dest='plot',
#                      help = 'Plot results.')
#    parser.add_option('-s', '--save', action='store_true', dest='save',
#                      help = 'Save the plots.')
    options, args = parser.parse_args()
    return options
            
def main():
    options = parse_args()
    data_00357 = load_data_stout('s', 0.00357, plist_stout, gflist_stout)
    data_01 = load_data_stout('s', 0.01, plist_stout, gflist_stout)
    
    #print [d.ap for d in data_00357]
#    print ar([d.bilinear_Lambda(1,0)[0].real for d in data_00357])
#    print ar([d.bilinear_Lambda(1,0)[1].real for d in data_00357])
#    print ''
#    print ar([d.bilinear_Lambda(2,0)[0].real for d in data_00357])
#    print ar([d.bilinear_Lambda(2,0)[1].real for d in data_00357])
#    print ''
#    print ar([d.bilinear_Lambda(4,0)[0].real for d in data_00357])
#    #print [d.bilinear_Lambda(4,0)[1].real for d in data_00357]
#    print ''
#    print ar([d.bilinear_Lambda(8,0)[0].real for d in data_00357])
#    #print [d.bilinear_Lambda(8,0)[1].real for d in data_00357]
#    print '\n\n'
    print ar([(2/d.Lambda_V(0)[0].real) for d in data_00357])
#    #print [d.Lambda_V(0)[1].real for d in data_00357]
    print ''
    print ar([(2/d.Lambda_V(1)[0].real) for d in data_00357])
#    #print [d.Lambda_V(1)[1].real for d in data_00357]
    print ''
    print ar([(2/d.Lambda_V(2)[0].real) for d in data_00357])
#    #print [d.Lambda_V(2)[1].real for d in data_00357]
    print ''
    print ar([(2/d.Lambda_V(3)[0].real) for d in data_00357])
    #print [d.Lambda_V(3)[1].real for d in data_00357]
    print ''
    print [d.Zq()[0].real for d in data_00357]
    print ''
    print [d.Zq2()[0].real for d in data_00357]
    
    print ''
    print [d.bilinear_Lambda(0,0)[0].real for d in data_00357]
    print [d.bilinear_Lambda(0,0)[1].real for d in data_00357]
    print [d.bilinear_Z(0,0)[0].real for d in data_00357]
    print [d.bilinear_Z(0,0)[1].real for d in data_00357]
    #print ''
    #print [d.Zq()[0].real for d in data_00357]
    #print [d.Zq2()[1].real for d in data_00357]
           
#    print [d.bilinear_Lambda(0,0)[0].real for d in data_01]
#    print [d.bilinear_Lambda(0,0)[1].real for d in data_01]
#    print ''
#    print [d.bilinear_Lambda(15,15)[0].real for d in data_01]
#    print [d.bilinear_Lambda(15,15)[1].real for d in data_01]

    #props_HYP = load_data_HYP('p', 0.01, plist_HYP, gflist_HYP)
    #props_naive = load_data_naive('n', 0.03, plist_naive, gflist_naive)
      
    #make_plots([data_00357, data_01], save=False)
    #compare_Zq([props_01, props_HYP, props_naive], save=False)
    
    return 0
    
def compare_Zq(data_list, save=False):
    '''Compare Zq using different actions.'''
    
    fontX = {'family':'monospace'}
    p.rc('font', **fontX)
    
    legend = ()
    p.figure()
    p.title('$\mathtt{Preliminary}$')
    p.xlabel('$\mathtt{\mu \quad [GeV]}$')
    p.ylabel('$\mathtt{Zq}$')
    for data in data_list:
        x = [d.mu for d in data]
        y, s = zip(*[d.Zq2() for d in data])
        dada = p.errorbar(x, y, s)
        legend += dada[0],
    p.legend(legend, ('stout', 'HYP', 'naive'), 'best')
    if save:
        root = '/Users/atlytle/Dropbox/TeX_docs/stout_NPR/figs/'
        p.savefig(root + 'Zq_compare.pdf')
    else:
        p.show()    

def plot_Zq(data_list, save=False):
    legend = ()
    p.figure()
    p.title(r'$\mathtt{Preliminary}$')
    p.xlabel('$\mathtt{[ap]^2}$')
    p.ylabel('$\mathtt{Zq}$')
    for data, f in zip(data_list, ['s','^']):
        x = [d.apSq for d in data]
        y, s = zip(*[d.Zq() for d in data])
        dada = p.errorbar(x, y, s, fmt='k'+f, mfc='none', ms=8, mew=1)
        legend += dada[0],
        y, s = zip(*[d.Zq2() for d in data])
        dada = p.errorbar(x, y, s, fmt='b'+f, mfc='none', mec='blue', ms=8, mew=1)
        legend += dada[0],
    p.legend(legend, ('$\mathtt{ap}$', '$\mathtt{sin[ap]}$'), 'best')
    if save:
        root = '/Users/atlytle/Dropbox/TeX_docs/stout_NPR/figs/'
        p.savefig(root + 'Zq_prelim.pdf')
    else:
        p.show()
        
def plot_M(data_list, save=False):
    legend = ()
    p.figure()
    p.title('$\mathtt{Preliminary}$')
    p.xlabel('$\mathtt{[ap]^2}$')
    p.ylabel('$\mathtt{M}$')
    for data, f in zip(data_list, ['ks', 'k^']):
        x = [d.apSq for d in data]
        y, s = zip(*[d.M() for d in data])
        dada = p.errorbar(x, y, s, fmt=f, mfc='none', ms=8, mew=1)
        legend += dada[0],
    p.legend(legend, ('am=0.00357', 'am=0.01'))
    if save:
        root = '/Users/atlytle/Dropbox/TeX_docs/stout_NPR/figs/'
        p.savefig(root + 'M_prelim.pdf')
    else:
        p.show()    
        
def plot_LambdaSP(data_list, save=False):
    legend = ()
    p.figure()
    p.title('$\mathtt{Preliminary}$')
    p.xlabel('$\mathtt{[ap]^2}$')
    p.ylabel('$\mathtt{\Lambda{S \otimes F}}$')
    for data, f in zip(data_list, ['s', '^']):
        x = [d.apSq for d in data]
        y, s = zip(*[d.bilinear_Lambda(0,0) for d in data])
        dada = p.errorbar(x, y, s, fmt='k'+f, mfc='none', ms=8, mew=1)
        legend += dada[0],
        y, s = zip(*[d.bilinear_Lambda(15,15) for d in data])
        dada = p.errorbar(x, y, s, fmt='r'+f, mec='red', mfc='none', ms=8, mew=1)
        legend += dada[0],
    p.legend(legend, ('$\mathtt{1 \otimes 1}$', 
                      '$\mathtt{\gamma_5 \otimes \gamma_5}$'))
    if save:
        root = '/Users/atlytle/Dropbox/TeX_docs/stout_NPR/figs/'
        p.savefig(root + 'LambdaSP_prelim.pdf')
    else:
        p.show()     
        
def make_plots(data_list, save=False):

    # Vintage look.
    fontX = {'family':'monospace'}
    p.rc('font', **fontX)

    # Plot Zq.
    plot_Zq(data_list, save)
    
    # Plot M.
    plot_M(data_list, save)
    
    # Plot Lambda_S,P.
    plot_LambdaSP(data_list, save)
    
if __name__ == "__main__":
    sys.exit(main())
