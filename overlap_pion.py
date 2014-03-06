import numpy as np

from read_overlap import correlator_name
from resample import yield_JK_sample_np, JK_block, block
from calc_overlap import plot_correlator, plot_effmass, fit_twopoint_cfuns

root = ''
mlist = [0.0745]

cs = [np.load(correlator_name(m)) for m in mlist]
csJK = map(JK_block, cs)
print cs[0].shape
print csJK[0].shape
t=0
for x in csJK[0][0]:  # correlator
    print t, x
    t += 1

#plot_correlator(csJK[0])
#plot_effmass(csJK[0])
print fit_twopoint_cfuns(csJK[0]/10000., 10, 22, 64)

