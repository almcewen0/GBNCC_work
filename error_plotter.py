import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rc("font",family='Arial')
fn={'fontname':'Arial'}
fs_lab=15
fs_tit=19
fs_leg=12

p0b=np.arange(28,200,1)
p0s=0.00600535*p0b-0.28950578
dc=0.08# np.arange(0.01,0.2,0.001)
w=dc*p0b
dw=1.#np.arange(0.1,15,0.1)
snm=10 # np.arange(10,100)
yoffsum=1.
sig=1. 
dsig=(snm/sig)*np.sqrt((p0s/120.)*yoffsum**2 + (dw/((2*sig)*(p0b-w)))**2)
dwer=snm*dw/w
dpi=(np.sqrt(p0s/120.))*p0b/np.sqrt(w)
dpo=p0b/(np.sqrt(w*(p0b-w)))
dsnm=np.sqrt(dpo**2+dsig**2+dwer**2+dpi**2)
plt.figure(figsize=(6.5,4.5))
xplt=p0b
plt.plot(xplt,(np.zeros(len(xplt))+dpi)/snm,label=r'$dp_{i}$',ls='--')
plt.plot(xplt,(np.zeros(len(xplt))+dpo)/snm,label=r'$dp_{off}$',ls='--')
plt.plot(xplt,(np.zeros(len(xplt))+dsig)/snm,label=r'$d\sigma_{off}$',ls='--')
plt.plot(xplt,(np.zeros(len(xplt))+dwer)/snm,label=r'$dW$',ls='--')
plt.plot(xplt,dsnm/snm,label=r'$dS/N$')
plt.legend(fontsize=fs_leg,loc='upper left')
plt.suptitle('Normalized error vs. Spin Period',fontsize=fs_tit)
plt.xlabel('Spin Period [bins]',fontsize=fs_lab)
plt.ylabel('dS/N / S/N',fontsize=fs_lab)
plt.show()
