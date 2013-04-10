import sys
import optparse
import read_staggered as rs
import staggered_NPR as sn

plist = [[2,2,2,4]]
gflist = [100, 200, 300, 400]

def load_data(ensemble, mass, plist, gflist):
    return [rs.StaggeredExceptional(ensemble, mass, p, gflist)
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
    
    
    return 0
    
if __name__ == "__main__":
    sys.exit(main())
