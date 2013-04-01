# Read-in and construct polespace propagators.

root = '/Users/atlytle/Documents/staggered_data/'

def prop_location(ensemble, mass, p, gf):
    '''Location of polespace propagators.'''
    ms = str(mass)
    ps = '.'.join(map(str, p))
    gfs = str(gf)
    loc = '{0}cluster_props/m{1}/{2}/data/pole_{2}_{3}.npy'.format(ensemble, 
                                                                   ms, ps, gfs)
    return root + loc
