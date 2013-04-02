import pylab as p
from parse_overlap import OverlapPoint, OverlapWall

def plot_correlator(cfnc):
    p.figure()
    p.xlabel('$t$')
    p.ylabel('$C[t]$')
    # for now just show central values
    p.semilogy(cfnc[0].real)
    p.show()
    
def plot_effmass(cfnc):
    pass
    
def naive_effmass(cfnc):
    '''Compute effective mass of correlation function.'''
    pass
    
    
def main():
    ms = 0.0495
    wall = OverlapWall(0.55, ms)
    point = OverlapPoint(0.551, ms)
    
    plot_correlator(wall.pscalar)

if __name__ == "__main__":
    main()
