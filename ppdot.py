import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from glob import glob

#figsize=(10,8)
#import matplotlib as mpl
#fs_lab=22
#fs_leg=16
#fs_tit=22
#fs_ticks=14
#mpl.rc('font',family='Arial')

import pylab

fontsize=14
fs_lab=20
fs_leg=14
pylab.rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
pylab.rc('xtick',labelsize=fontsize)
pylab.rc('ytick',labelsize=fontsize)
pylab.rc('text',usetex=True)


res=glob('master_results_*.txt')
if len(res)>1:
  diff=[]
  import datetime as dt
  for f in res:
      tim=str(f.split('.')[0].split('_')[-1])
      dat=dt.date(month=int(tim[0:2]),day=int(tim[2:4]),year=int(tim[4:]))
      diff.append((dt.date.today()-dat).days)
  diff=np.array(diff)
  res=np.array(res)[diff==diff.min()][0]
elif len(res)==1:
  res=res[0]
else:
  print "missing results file"
  exit()
resn=Table.read(res,format='ascii')
full=Table.read('/Users/aemcewen/GBNCC_paper/psr_lists/full_list.txt',format='ascii')

cond=[]

plt.figure()#figsize=figsize)
plt.yscale('log')
plt.xscale('log')
ax=plt.axes()
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
condnd=[p['name'] not in resn['name'] and p['p1']!='*' for p in full]
x=np.arange(5e-4,12,1)
tau=np.array([10e3,100e3,1e6,10e6,100e6,1e9,10e9,100e9])
for tc in tau:
    y=15.8e-9*x/tc
    plt.plot(x,y,c='black',alpha=0.5,ls='--',lw=1)
    tct=int(("%1.0e" %tc).split('+')[-1])
    plt.text(x[-1],2.5*y[-5],'10$^{%s}$ yr' %tct,rotation=13,bbox=dict(facecolor='white',edgecolor='white'),fontsize=fontsize)
bc=np.array([1e8,1e9,1e10,1e11,1e12,1e13,1e14])
i=0
x2=[2.5e-2,0.2,2.5,7e-4,5e-3,5e-2,0.6]
y2=[1e-22,1.25e-21,1e-20,5e-15,4.5e-14,5e-13,4.55e-12]
for b in bc:
    y=10**(-15)*(b/1e12)**2/x
    plt.plot(x,y,c='black',alpha=0.5,ls='-.',lw=1)
    tct=int(("%1.0e" %b).split('+')[-1]) 
    plt.text(x2[i],2.5*y2[i],'10$^{%s}$ G' %tct,rotation=-15,bbox=dict(facecolor='white',edgecolor='white'),fontsize=fontsize)
    i+=1
plt.plot(full[condnd]['p0'],full[condnd]['p1'].astype(float),marker='o',markersize=6,alpha=0.5,color='grey',linewidth=0,label='Non-detections')
plt.plot(resn[resn['p1']!='*']['p0'],resn[resn['p1']!='*']['p1'].astype(float),marker='+',markersize=9,linewidth=0,label='Detections')
plt.xlim([6e-4,11])
plt.ylim([5e-23,1e-10])
plt.legend(loc='upper left',fontsize=fs_leg)
plt.xlabel('Spin Period [s]',fontsize=fs_lab)
plt.ylabel('Period Derivative [s s$^{-1}$]',fontsize=fs_lab)
plt.tick_params()#labelsize=fs_ticks)
plt.subplots_adjust(left=0.14,bottom=0.12,right=0.9,top=0.98)
plt.savefig('plots/ppdot.pdf')
#plt.show()


