import sys
import subprocess
import os
import numpy as np
from xml.dom import minidom

def parse_xml(elem):
    '''Convert each elem to a complex number.'''
    re = float(elem.getElementsByTagName('re')[0].firstChild.data)
    im = float(elem.getElementsByTagName('im')[0].firstChild.data)
    return complex(re,im)
    
def create_npy(dat):
    '''Create polespace propagator data structure.'''
    
    ar = np.array(dat)  # One dimensional.
    print ar.shape
    ar = ar.reshape(16,16,3,3)  # A_{ijab}; rightmost index "fastest".
    ar = ar.transpose((0,2,1,3)) # A_{iajb}.
    ar = ar.reshape(48,48) # A_{ia,jb} (Direct product form.)
    
    return ar
    

def main(argv):
    '''Parse polespace xml into a numpy array.'''
    
    # Create data/ subdirectory if it doesn't already exist.
    if not os.path.isdir('./data/'):
        print 'Creating data/ directory...'
        subprocess.call('mkdir data', shell=True)

    for arg in argv:
        basestring = arg.replace('.xml', '')
        target = './data/{0}.npy'.format(basestring)
        if os.path.isfile(target):
            print target, 'already exists!'
            continue

        xmldoc = minidom.parse(arg) #xml check?
        cnums = map(parse_xml, xmldoc.getElementsByTagName('elem'))

        a = create_npy(cnums)
        print a

        print arg, 'has been parsed.'
        
    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
