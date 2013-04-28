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

        self.prop_list = [np.load(self.prop_location(gf))
                          for gf in gflist]
        
        self.prop, self.propJK = npr.JKcompute(lambda x: np.average(x, axis=0),
                                               self.prop_list)
                                               
        # Initialize bilinear containers.
        self.bilinear_list = np.array([[None]*16]*16)
        self.bilinear = np.array([[None]*16]*16)
        self.bilinearJK = np.array([[None]*16]*16)
        
        self.bilinear_list[0,0] = [np.load(self.bilinear_location(0,0,gf))
                                   for gf in gflist]
        self.bilinear[0][0], self.bilinearJK[0][0] = npr.JKcompute(
                lambda x: np.average(x, axis=0), self.bilinear_list[0][0])
        self.bilinear_list[15,15] = [np.load(self.bilinear_location(15,15,gf))
                                     for gf in gflist]
        self.bilinear[15][15], self.bilinearJK[15][15] = npr.JKcompute(
                lambda x: np.average(x, axis=0), self.bilinear_list[15][15])


    def prop_location(self, gf):
        '''Location of polespace propagators.'''
        name = 'pole_{0}_{1}.npy'.format(self.pstring, gf)
        
        return self.root + name
        
    def bilinear_location(self, S, F, gf):
        '''Location of bilinear correlation function SxF.'''
        name = 'bl_{0}_{1}_{2}_{3}.npy'.format(self.pstring, gf, S, F)
        
        return self.root + name
        
        
    def populate_kinematic_variables(self):
        '''Assign momentum variables: ap, mu, etc.'''
        p, L, T, a = self.p, self.L, self.T, self.a
        self.ap = (2*np.pi)*np.array([p[0]/L, p[1]/L, p[2]/L,
                                      p[3]/T])
        self.apSq = np.dot(self.ap, self.ap)
        self.mu = np.sqrt(self.apSq)/a
        self.apHat = map(np.sin, self.ap)
                
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
        
    def Zq2(self):
        '''Zq calculated from the propagator (RI' scheme).'''
        ap = np.sin(self.ap)
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
    
    def bilinear_Lambda(self, S, F):
        '''Trace of amputated bilinear correlation function.'''
        calc = lambda prop, bl :npr.bilinear_Lambda(prop, bl, S, F)/self.V
        r, JK = calc(self.prop, self.bilinear[S][F]),\
                [calc(*x) for x in zip(self.propJK, self.bilinearJK[S][F])]
        return r, npr.JKsigma(r, JK)
        
class StoutExceptional(StaggeredExceptional):
    a = 1/1.09  # 1/GeV.
    L, T = 24., 32.
    V = (L**3)*T
    def __init__(self, ensemble, mass, p, gflist):
        StaggeredExceptional.__init__(self, ensemble, mass, p, gflist)
        self.populate_kinematic_variables()
        
class HYPExceptional(StaggeredExceptional):
    a = 0.602  # 1/GeV. (corresponds to am_sea = .01 ensemble)
    L, T = 20., 64.
    V = (L**3)*T
    def __init__(self, ensemble, mass, p, gflist):
        StaggeredExceptional.__init__(self, ensemble, mass, p, gflist)
        self.populate_kinematic_variables() 
