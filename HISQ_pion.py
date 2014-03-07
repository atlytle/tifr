import numpy as np

from read_HISQ import correlator_name
from resample import JK_block

root = ''
mlist = ['635', '0102', '0509']
cs = [np.load(correlator_name(m)) for m in mlist]
csJK = map(JK_block, cs)

print cs[0].shape

t = 0
for x in cs[0][0]:
    print t, x
    t += 1
