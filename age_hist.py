import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt
from glob import glob
import pylab

fontsize=14
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

p1=np.array(resn[resn['p1']!='*']['p1'].astype(float)) 
p0=np.array(resn[resn['p1']!='*']['p0'].astype(float))
ages=np.array(15.8e-9*p0/p1)
log=raw_input("log/lin (l/i)?")
if log=='l':
    log=True
else:
    log=False

save=raw_input("save (y/n)?")
if save=='y':
    save=True
else:
    save=False
fs=10
plt.clf()

if log:
    plt.yscale('log')
    pos=np.logspace(np.log10(24),np.log10(40),4)
else:
    pos=np.linspace(33,40,4)

n,b,p=plt.hist(np.log10(ages),bins=30)
plt.text(9,pos[3],"N$_{psr}$ = %i" %len(ages),fontsize=fs)
plt.text(9.3,pos[2],"$\mu$ = %2.3g yr" %10**b.mean(),fontsize=fs)
plt.text(9.3,pos[1],"$\sigma$ = %2.3g yr" %10**b.std(),fontsize=fs)
plt.text(8.75,pos[0],"median = %2.3g yr" %10**np.median(b),fontsize=fs)
plt.axvline(x=b.mean(),c='red')
plt.axvline(x=b.mean()-b.std()/2,ls='--',c='red')
plt.axvline(x=b.mean()+b.std()/2,ls='--',c='red')
plt.xlabel('Age [log(yr)]',fontsize=12)
plt.ylabel('Count',fontsize=12)

if save:
    if log:
        plt.savefig('plots/gbncc_agehist_log.png')
    else:
        plt.savefig('plots/gbncc_agehist_linear.png')
else:
    plt.show()
