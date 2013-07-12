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
        with open(self.loc + 'a4a4_sc.dat') as f:
            self.a4a4 = jackknife_correlators(parse_correlator(f))
        with open(self.loc + 'psa4_sc.dat') as f:
            self.psa4 = jackknife_correlators(parse_correlator(f))
        with open(self.loc + 'a4ps_sc.dat') as f:
            self.a4ps = jackknife_correlators(parse_correlator(f))
            
        assert self.pscalar.shape == self.vector.shape == self.a4a4.shape == \
               self.psa4.shape == self.a4ps.shape
        self.N, self.T = self.pscalar.shape
            
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
        with open(self.loc + 'a4a4_sc.dat') as f:
            self.a4a4 = jackknife_correlators(parse_correlator(f))
        with open(self.loc + 'psa4_sc.dat') as f:
            self.psa4 = jackknife_correlators(parse_correlator(f))
        with open(self.loc + 'a4ps_sc.dat') as f:
            self.a4ps = jackknife_correlators(parse_correlator(f))
            
        assert self.pscalar.shape == self.vector.shape == self.a4a4.shape == \
               self.psa4.shape == self.a4ps.shape
        self.N, self.T = self.pscalar.shape

class Overlap_T96:
    def __init__(self, mc, ms):
        self.mc, self.ms = mc, ms
        self.root = '/Users/atlytle/Documents/overlap/L32T96/'
        self.load()
        
    def load(self):
        # Load point-point and point-wall correlators.
        self.loc = self.root + 'point/COR_c{0}_s{1}/'.format(self.mc, self.ms)
        
        with open(self.loc + 'pscalar_sc.dat') as f:
            self.pscalar_pp = jackknife_correlators(parse_correlator(f))
        with open(self.loc + 'vector_sc.dat') as f:
            self.vector_pp = jackknife_correlators(parse_correlator(f))
        with open(self.loc + 'a4a4_sc.dat') as f:
            self.a4a4_pp = jackknife_correlators(parse_correlator(f))
        with open(self.loc + 'psa4_sc.dat') as f:
            self.psa4_pp = jackknife_correlators(parse_correlator(f))
        with open(self.loc + 'a4ps_sc.dat') as f:
            self.a4ps_pp = jackknife_correlators(parse_correlator(f))
            
        assert self.pscalar_pp.shape == self.vector_pp.shape == \
               self.a4a4_pp.shape == self.psa4_pp.shape == self.a4ps_pp.shape
        self.N, self.T = self.pscalar_pp.shape
        
        with open(self.loc + 'p2w_pscalar_sc.dat') as f:
            self.pscalar_pw = jackknife_correlators(parse_correlator(f))
        with open(self.loc + 'p2w_vector_sc.dat') as f:
            self.vector_pw = jackknife_correlators(parse_correlator(f))
        with open(self.loc + 'p2w_a4a4_sc.dat') as f:
            self.a4a4_pw = jackknife_correlators(parse_correlator(f))
        with open(self.loc + 'p2w_psa4_sc.dat') as f:
            self.psa4_pw = jackknife_correlators(parse_correlator(f))
        with open(self.loc + 'p2w_a4ps_sc.dat') as f:
            self.a4ps_pw = jackknife_correlators(parse_correlator(f))
        
        assert self.pscalar_pw.shape == self.vector_pw.shape == \
               self.a4a4_pw.shape == self.psa4_pw.shape == self.a4ps_pw.shape
        assert self.N, self.T == self.pscalar_pw.shape
        
        # Load wall-wall and wall-point correlators.
        self.loc = self.root + 'wall/COR_c{0}_s{1}/'.format(self.mc, self.ms)
        
        with open(self.loc + 'pscalar_sc.dat') as f:
            self.pscalar_wp = jackknife_correlators(parse_correlator(f))
        with open(self.loc + 'vector_sc.dat') as f:
            self.vector_wp = jackknife_correlators(parse_correlator(f))
        with open(self.loc + 'a4a4_sc.dat') as f:
            self.a4a4_wp = jackknife_correlators(parse_correlator(f))
        with open(self.loc + 'psa4_sc.dat') as f:
            self.psa4_wp = jackknife_correlators(parse_correlator(f))
        with open(self.loc + 'a4ps_sc.dat') as f:
            self.a4ps_wp = jackknife_correlators(parse_correlator(f))
            
        assert self.pscalar_wp.shape == self.vector_wp.shape == \
               self.a4a4_wp.shape == self.psa4_wp.shape == self.a4ps_wp.shape
        assert self.N, self.T == self.pscalar_wp.shape
        
        with open(self.loc + 'p2w_pscalar_sc.dat') as f:
            self.pscalar_ww = jackknife_correlators(parse_correlator(f))
        with open(self.loc + 'p2w_vector_sc.dat') as f:
            self.vector_ww = jackknife_correlators(parse_correlator(f))
        with open(self.loc + 'p2w_a4a4_sc.dat') as f:
            self.a4a4_ww = jackknife_correlators(parse_correlator(f))
        with open(self.loc + 'p2w_psa4_sc.dat') as f:
            self.psa4_ww = jackknife_correlators(parse_correlator(f))
        with open(self.loc + 'p2w_a4ps_sc.dat') as f:
            self.a4ps_ww = jackknife_correlators(parse_correlator(f))
        
        assert self.pscalar_ww.shape == self.vector_ww.shape == \
               self.a4a4_ww.shape == self.psa4_ww.shape == self.a4ps_ww.shape
        assert self.N, self.T == self.pscalar_ww.shape
        
class Overlap_T144:
    def __init__(self, mc, ms):
        self.mc, self.ms = mc, ms
        self.root = '/Users/atlytle/Documents/overlap/L48T144/'
        self.load()
        
    def load(self):
        # Load point-point and point-wall correlators.
        self.loc = self.root + 'point/COR_c{0}_s{1}/'.format(self.mc, self.ms)
        
        with open(self.loc + 'pscalar_sc.dat') as f:
            self.pscalar_pp = jackknife_correlators(parse_correlator(f))
        with open(self.loc + 'vector_sc.dat') as f:
            self.vector_pp = jackknife_correlators(parse_correlator(f))
        with open(self.loc + 'a4a4_sc.dat') as f:
            self.a4a4_pp = jackknife_correlators(parse_correlator(f))
        with open(self.loc + 'psa4_sc.dat') as f:
            self.psa4_pp = jackknife_correlators(parse_correlator(f))
        with open(self.loc + 'a4ps_sc.dat') as f:
            self.a4ps_pp = jackknife_correlators(parse_correlator(f))
            
        assert self.pscalar_pp.shape == self.vector_pp.shape == \
               self.a4a4_pp.shape == self.psa4_pp.shape == self.a4ps_pp.shape
        self.N, self.T = self.pscalar_pp.shape
        
        with open(self.loc + 'p2w_pscalar_sc.dat') as f:
            self.pscalar_pw = jackknife_correlators(parse_correlator(f))
        with open(self.loc + 'p2w_vector_sc.dat') as f:
            self.vector_pw = jackknife_correlators(parse_correlator(f))
        with open(self.loc + 'p2w_a4a4_sc.dat') as f:
            self.a4a4_pw = jackknife_correlators(parse_correlator(f))
        with open(self.loc + 'p2w_psa4_sc.dat') as f:
            self.psa4_pw = jackknife_correlators(parse_correlator(f))
        with open(self.loc + 'p2w_a4ps_sc.dat') as f:
            self.a4ps_pw = jackknife_correlators(parse_correlator(f))
        
        assert self.pscalar_pw.shape == self.vector_pw.shape == \
               self.a4a4_pw.shape == self.psa4_pw.shape == self.a4ps_pw.shape
        assert self.N, self.T == self.pscalar_pw.shape
        
        # Load wall-wall and wall-point correlators.
        self.loc = self.root + 'wall/COR_c{0}_s{1}/'.format(self.mc, self.ms)
        
        with open(self.loc + 'pscalar_sc.dat') as f:
            self.pscalar_wp = jackknife_correlators(parse_correlator(f))
        with open(self.loc + 'vector_sc.dat') as f:
            self.vector_wp = jackknife_correlators(parse_correlator(f))
        with open(self.loc + 'a4a4_sc.dat') as f:
            self.a4a4_wp = jackknife_correlators(parse_correlator(f))
        with open(self.loc + 'psa4_sc.dat') as f:
            self.psa4_wp = jackknife_correlators(parse_correlator(f))
        with open(self.loc + 'a4ps_sc.dat') as f:
            self.a4ps_wp = jackknife_correlators(parse_correlator(f))
            
        assert self.pscalar_wp.shape == self.vector_wp.shape == \
               self.a4a4_wp.shape == self.psa4_wp.shape == self.a4ps_wp.shape
        assert self.N, self.T == self.pscalar_wp.shape
        
        with open(self.loc + 'p2w_pscalar_sc.dat') as f:
            self.pscalar_ww = jackknife_correlators(parse_correlator(f))
        with open(self.loc + 'p2w_vector_sc.dat') as f:
            self.vector_ww = jackknife_correlators(parse_correlator(f))
        with open(self.loc + 'p2w_a4a4_sc.dat') as f:
            self.a4a4_ww = jackknife_correlators(parse_correlator(f))
        with open(self.loc + 'p2w_psa4_sc.dat') as f:
            self.psa4_ww = jackknife_correlators(parse_correlator(f))
        with open(self.loc + 'p2w_a4ps_sc.dat') as f:
            self.a4ps_ww = jackknife_correlators(parse_correlator(f))
        
        assert self.pscalar_ww.shape == self.vector_ww.shape == \
               self.a4a4_ww.shape == self.psa4_ww.shape == self.a4ps_ww.shape
        assert self.N, self.T == self.pscalar_ww.shape
        


