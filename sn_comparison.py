import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from astropy.table import Table
from glob import glob
import datetime as dt
import matplotlib as mpl

#figsize=(10,8)
#import matplotlib as mpl
#fs_lab=22
#fs_leg=16
#fs_tit=22
#fs_ticks=16
#mpl.rc('font',family='Arial')

import pylab

fontsize=14
fs_lab=20
pylab.rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
pylab.rc('xtick',labelsize=fontsize)
pylab.rc('ytick',labelsize=fontsize)
pylab.rc('text',usetex=True)

def dot(l1,l2):
      l1=l1.astype(int)
      l2=l2.astype(int)
      if len(l1)==len(l2):
          l=l1+l2-1
          for i in range(len(l)):
             if l[i]!=1:
                l[i]=0
      return l.astype(bool)

def snm_err(x,y,p0s,snm,w,dw,lim):
    p0b=len(x)
    dsig=(snm/y[lim].std())*np.sqrt((p0s/120.)*np.sum(y[lim]-y[lim].mean())**2 + (dw/(2*y[lim].std()*(p0b-w)))**2)
    dwer=snm*dw/w
    dpi=(np.sqrt(p0s/120.))*p0b/np.sqrt(w)
    dpo=p0b/(np.sqrt(w*(p0b-w)))
    dsnm=np.sqrt(dpo**2+dsig**2+dwer**2+dpi**2)
    return dsnm



res=glob('master_results_*.txt')
if len(res)>1:
  import datetime as dt
  diff=[]
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
cond=resn['sne']!='*' 


cm=plt.cm.get_cmap('RdYlBu') 
fig,ax=plt.subplots()#figsize=figsize)
ax.set_xlabel('Expected S/N',fontsize=fs_lab)
ax.set_ylabel('Measured S/N',fontsize=fs_lab)
ax.set_xscale('log')
ax.set_yscale('log')
#fig.suptitle('Measured vs. Expected S/N in GBNCC')#,fontsize=fs_tit)
colorax='p0'
try: 
    len(resn['%s' %colorax]!='*')>1
    cond=dot(cond,resn['%s' %colorax]!='*')
except:
    pass
x=resn[cond]['sne'].astype(float)
y=resn[cond]['snm']
z=resn[cond]['%s' %colorax].astype(float)
dy=[]
for l in resn[cond]:
    nm=l['name']
    p0s=l['p0']
    w=l['on-pulse']
    snm=l['snm']
    try:
        prof=glob('./presto_plots/master_112818/good_profs_withwidth/%s*bestprof*' %nm)[0]
    except:
	print "no prof for "+nm
	continue
    n,v,lim=np.loadtxt(prof,unpack=True)
    lim=lim.astype(bool)
    dy.append(snm_err(n,v,p0s,snm,w,1.0,lim))
#sc=ax.scatter(x,y,c=np.log10(z),cmap=cm,marker='+')
ln=[2,8e2]
ax.set_aspect('equal')
ax.plot(ln,ln)
ax.set_xscale('log')
ax.set_yscale('log')
ax.plot(resn[cond]['sne'].astype(float),resn[cond]['snm'],c='blue',marker='+',lw=0)
ax.errorbar(resn[cond]['sne'].astype(float),resn[cond]['snm'],yerr=dy/resn['snm'][cond],lw=0,elinewidth=1,capsize=2,c='blue')

#cb=plt.colorbar(sc,ticks=[])
#cb.set_label('Spin Period',fontsize=fs_lab)
plt.subplots_adjust(left=0.10,right=0.98,top=0.92,bottom=0.13,wspace=0.0,hspace=0.0)
dat=str((dt.datetime.now().date())).split('-') 
curdat=dat[1]+dat[2]+dat[0][2:]; 
plt.tick_params()#labelsize=fs_ticks)
plt.savefig('plots/snm_vs_sne_%s.pdf' %curdat)
#plt.show()
