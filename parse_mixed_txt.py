"""Parse pion data out of the correlator text files from Subhasis.

Returns the pion data as text files, as opposed to .npy format
in parse_mixed.py.  This makes it easier to do things like correlated
jackknife fits.
"""

import sys

def in_name(m1, m2, config):
    """Input text file."""
    root = '/Users/atlytle/Documents/dmix/for_tifr/corr'
    return root + "mix_ksm{0}_ovm{1}_corr.{2}".format(m1, m2, config)

def main():
    pass

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))