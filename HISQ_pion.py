import numpy as np

from read_HISQ import correlator_name
from resample import JK_block
from calc_overlap import plot_correlator, plot_effmass, fit_twopoint_cfuns

root = ''
mlist = ['635', '0102', '0509']
cs = [np.load(correlator_name(m)) for m in mlist]
csJK = map(JK_block, cs)

print cs[0].shape

t = 0
for x in csJK[0][0]:
    print t, x
    t += 1
    
#plot_correlator(csJK[0])
#plot_correlator(csJK[1])
#plot_correlator(csJK[2])

#plot_effmass(csJK[0])
#plot_effmass(csJK[1])
#plot_effmass(csJK[2])

print fit_twopoint_cfuns(csJK[0]/10000., 10, 22, 64)
print fit_twopoint_cfuns(csJK[1]/100000., 10, 22, 64)
print fit_twopoint_cfuns(csJK[2]/10000., 10, 22, 64)
