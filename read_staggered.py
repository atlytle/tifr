import numpy as np

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

    def prop_location(self, gf):
        '''Location of polespace propagators.'''
        name = 'pole_{0}_{1}.npy'.format(self.pstring, gf)
        
        return self.root + name
    
