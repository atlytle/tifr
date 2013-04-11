import sys
import optparse
import pylab as p
import StaggeredExceptional as se
import staggered_NPR as sn

plist = [(1,1,1,2), (2,2,2,2), (2,2,2,3), (2,2,2,4), (3,2,2,3),
         (2,3,3,3), (3,2,3,4), (3,3,3,3)]
gflist = [100, 200, 300, 400]

def load_data(ensemble, mass, plist, gflist):
    return [se.StaggeredExceptional(ensemble, mass, p, gflist)
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
    props = load_data('s', 0.00357, plist, gflist)
    #print [p.M() for p in props]
    #print ''
    #print [p.Zq()[0] for p in props]
    
    make_plots(props)
    
    return 0
    
def make_plots(Data):
    p.figure()
    p.title(r'$\mathtt{Preliminary}$')
    p.xlabel('$\mathtt{(ap)^2}$')
    p.ylabel('$\mathtt{Z_q}$')
    x = [d.apSq for d in Data]
    y, s = zip(*[d.Zq() for d in Data])
    p.errorbar(x, y, s, fmt='ks', mfc='none', ms=8)
    p.show()
    
    p.figure()
    p.title('$\mathtt{Preliminary}$')
    p.xlabel('$(ap)^2$')
    p.ylabel('$\mathtt{M}$')
    x = [d.apSq for d in Data]
    y, s = zip(*[d.M() for d in Data])
    p.errorbar(x, y, s, fmt='k^', mfc='none', ms=8)
    p.show()
    
    
if __name__ == "__main__":
    sys.exit(main())
