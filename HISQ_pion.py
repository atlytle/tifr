import numpy as np

from read_HISQ import correlator_name, correlator_name2
from resample import JK_block
from calc_overlap import plot_correlator, plot_effmass, fit_twopoint_cfuns
from correlators import fit_cfuns_double_cosh_osc

#root = ''
#mlist = ['635', '0102', '0509']
#cs = [np.load(correlator_name(m)) for m in mlist]
#csJK = map(JK_block, cs)

#c_635_0509 = np.load(correlator_name2('635', '0509'))
#c_635_0509JK = JK_block(c_635_0509)


c_0509_635 = np.load(correlator_name2('0509', '635'))
c_0509_635JK = JK_block(c_0509_635)

c_0509_0509 = np.load(correlator_name2('0509', '0509'))
c_0509_0509JK = JK_block(c_0509_0509)


print c_0509_0509.shape

#t = 0
#for x in c_0509_635JK[0]:
#    print t, x.real,' ', x.imag
#    t += 1
    
    
plot_correlator(c_0509_635JK.real)
plot_effmass(c_0509_635JK.real)
print fit_twopoint_cfuns(c_0509_0509JK.real/10000., 10, 22, 64)
print fit_cfuns_double_cosh_osc(c_0509_0509JK.real/10000., 10, 22, 64)



#plot_correlator(c_0509_635JK.real)
#plot_effmass(c_0509_635JK.real)
#print fit_twopoint_cfuns(c_0509_635JK.real/10000., 10, 22, 64)

#plot_correlator(csJK[0])
#plot_correlator(csJK[1])
#plot_correlator(csJK[2])

#plot_effmass(csJK[0])
#plot_effmass(csJK[1])
#plot_effmass(csJK[2])

#print fit_twopoint_cfuns(csJK[0]/10000., 10, 22, 64)
#print fit_twopoint_cfuns(csJK[1]/100000., 10, 22, 64)
#print fit_twopoint_cfuns(csJK[2]/10000., 10, 22, 64)
