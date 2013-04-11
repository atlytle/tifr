import numpy as np
import staggered_NPR as npr
from numpy.linalg import inv

root = '/Users/atlytle/Documents/staggered_data/'

class StaggeredExceptional:
    def __init__(self, ensemble, mass, p, gflist):
        self.ensemble, self.m, self.p, self.gflist = ensemble, mass, p, gflist
        self.pstring = '.'.join(map(str, p))
        self.root = '/Users/atlytle/Documents/staggered_data/'\
                    '{0}cluster_props/m{1}/{2}/data/'.format(ensemble, str(mass),
                                                             self.pstring)
                                                             
        self.L, self.T =24., 32.
        self.ap = (2*np.pi)*np.array([p[0]/self.L, p[1]/self.L, p[2]/self.L,
                                      p[3]/self.T])
        self.apSq = np.dot(self.ap, self.ap)
        self.apHat = map(np.sin, self.ap)
        self.prop_list = [np.load(self.prop_location(gf))
                          for gf in gflist]
        
        self.prop, self.propJK = npr.JKcompute(lambda x: np.average(x, axis=0),
                                               self.prop_list)

    def prop_location(self, gf):
        '''Location of polespace propagators.'''
        name = 'pole_{0}_{1}.npy'.format(self.pstring, gf)
        
        return self.root + name
        
    def Zq(self):
        '''Zq calculated from the propagator (RI' scheme).'''
        ap = self.ap
        def calc_Zq(prop):
            Sinv = (1./2)*inv(prop)  # 1/2 from action def.
            tr = npr.ps_trace
            traces = tr(1,0,Sinv), tr(2,0,Sinv), tr(4,0,Sinv), tr(8,0,Sinv)
            Zinv = np.dot(ap, traces).imag/np.dot(ap, ap)
            return 1./Zinv
        r, JK = calc_Zq(self.prop), map(calc_Zq, self.propJK)
        return r, npr.JKsigma(r, JK)
        
        
        
    def M(self):
        '''Mass term calculated from the propagator.'''
        def calc_M(prop):
            Sinv = (1./2)*inv(prop)  # 1/2 from action def.
            return npr.ps_trace(0,0,Sinv).real  
        r, JK = calc_M(self.prop), map(calc_M, self.propJK)
        return r, npr.JKsigma(r, JK)
    
