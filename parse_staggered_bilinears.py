import sys
import subprocess
import os
import xml.etree.ElementTree as ET
import numpy as np

from multiprocessing import Pool

def parse_elem(elem):
    "Convert each elem to a complex number."
    re = float(elem.find('re').text)
    im = float(elem.find('im').text)
    return complex(re,im)
    
def create_npy(dat):
    '''Create bilinear data structure.'''
    
    ar = np.array(dat)  # One dimensional.
    ar = ar.reshape(16,16,3,3)  # A_{ijab}; rightmost index "fastest".
    ar = ar.transpose((0,2,1,3)) # A_{iajb}.
    ar = ar.reshape(48,48) # A_{ia,jb} (Direct product form.)
    
    return ar

def parse_xml(xml):
    "Memory efficient parsing routine extracts bilinear data from xml."

    print 'Parsing {0}...'.format(xml)
    basestring = xml.replace('.xml','')
    suppress, sup_limit = 0, 4
    for event, node in ET.iterparse(xml):

        tmp = node.tag.partition('_')[0]
        if tmp == 'bl' or tmp == 'blNE':  # Look for bl_S_F tags.
            tmp = node.tag.partition('_')[0]
            tag = node.tag.lstrip(tmp)  # _S_F identifier.
            target = './data/{0}{1}'.format(basestring, tag)
            if os.path.isfile(target+'.npy'):
                suppress += 1
                if suppress < sup_limit:
                    print target, 'already exists.'
                if suppress == sup_limit:
                    print target, 'already exists...suppressing output.'
                node.clear()  # Release node from memory.
                continue
            
            cnums = map(parse_elem, node.getiterator('elem'))
            node.clear()  # Release node from memory.
         
            a = create_npy(cnums)
            np.save(target, a)
    
    return 0

def main(argv):
    if not os.path.isdir('data'):
        os.mkdir('data')

    pool = Pool()
    for arg in argv:
        pool.apply_async(parse_xml, (arg,))
    pool.close()
    pool.join()
    
    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
