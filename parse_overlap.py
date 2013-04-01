import numpy as np

def parse_correlator(data):
    "Convert text correlator file into a numpy array."
    
    result = []
    time = []  # Temp variable to infer Nt.
    
    # Strip out numbers.
    for line in data:
        tmp = line.strip().split(' ')
        time.append(int(tmp[-3]))
        re, im = float(tmp[-2]), float(tmp[-1])
        z = complex(re, im)
        result.append(z)
    
    L = len(result)
    Nt = max(time) + 1
    assert L % Nt == 0  # Make sure Nconf is an integer.
    result = np.array(result)
    result = np.resize(result, (L/Nt, Nt))  # Nconf x Nt.
    return result
    
def jackknife_correlators(corrs):
    "Return average and jackknife averages of correlators."
    test = corrs.copy()
    Nconf = corrs.shape[0]  # Number of configs.
    jackknifes = corrs.sum(axis=0)*(1./Nconf)  # First row is average.
    for dummy in range(Nconf):
        c, cs = corrs[0:1], corrs[1:]
        jackknifes = np.vstack((jackknifes, cs.sum(axis=0)*(1./(Nconf-1))))
        corrs = np.vstack((cs,c))
    assert (test == corrs).all()
    return jackknifes
    
class OverlapWall:
    
    def __init__(self, mc, ms):
        self.mc, self.ms = mc, ms
        self.root = '/Users/atlytle/Documents/overlap/wall/'
        self.loc = self.root + 'COR_c{0}_s{1}/'.format(mc, ms)
        self.load()
        
    def load(self):
        with open(self.loc + 'pscalar_sc.dat') as f:
            self.pscalar = jackknife_correlators(parse_correlator(f))
        with open(self.loc + 'vector_sc.dat') as f:
            self.vector = jackknife_correlators(parse_correlator(f))
            
class OverlapPoint:
    
    def __init__(self, mc, ms):
        self.mc, self.ms = mc, ms
        self.root = '/Users/atlytle/Documents/overlap/point/'
        self.loc = self.root + 'COR_c{0}_s{1}/'.format(mc, ms)
        self.load()
        
    def load(self):
        with open(self.loc + 'pscalar_sc.dat') as f:
            self.pscalar = jackknife_correlators(parse_correlator(f))
        with open(self.loc + 'vector_sc.dat') as f:
            self.vector = jackknife_correlators(parse_correlator(f))
        


