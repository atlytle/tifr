import optparse
import sys
import itertools

from parse_overlap import Overlap_T96, Overlap_T144
from calc_overlap import plot_correlator, plot_effmass

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
    
    # L32T96.
    charm_masses = [0.425, 0.430, 0.435]
    strange_masses = [0.046, 0.048, 0.049]
    
    for mc, ms in itertools.product(charm_masses, strange_masses):
        print mc, ms
        x = Overlap_T96(mc, ms)
        name = root + 'L32T96_c{0}_u{1}_pscalar_pp_'.format(mc, ms)
        plot_correlator(x.pscalar_pp, options.save, name+'corr.pdf')
        plot_effmass(x.pscalar_pp, options.save, name+'meff_naive.pdf')
    
    return 0
    
if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
