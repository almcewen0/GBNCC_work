import numpy as np
import matplotlib.pyplot as plt
from glob import glob
from astropy.table import Table
from subprocess import Popen,PIPE
from functions_python import shift

#figsize=(10,8)
#import matplotlib as mpl
#fs_lab=20
#fs_leg=14
#fs_tit=20
#fs_ticks=12
fs_txt=7
#mpl.rc('font',family='Arial')

import pylab
fontsize=14
fs_lab=18
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
profs=[]
for p in resn['name']:
    ls=glob('presto_plots/master_112818/good_profs_withwidth/%s*bestprof*' %p)
    if len(ls)>1:
        print "Multiple plots for %s" %p
    profs.append(glob('presto_plots/master_112818/good_profs_withwidth/%s*bestprof*' %p)[0])
down=6
across=10
cent_freq=0.350 #GHz
df=0.0244140625 #MHz

I=len(profs)/(across*down)
J=across
K=down
for i in range(I):
    fig,ax=plt.subplots(across,down,figsize=(6,8))
    for j in range(J):
        for k in range(K):
            prof_ind=k+j*K+i*(J*K)
            n,v,lim=np.loadtxt(profs[prof_ind],unpack=True)
            v=(v/v.max()).tolist()
	    while v.index(np.array(v).max())!=int(n[-1]*3/4):
		v=shift(v,1)
	    v=np.array(v)
            nm=profs[prof_ind].split('/')[-1].split('_')[0]
            dm=resn[resn['name']==nm]['dms']
	    name=nm
	    if len(name.split('-'))>1:
		name=name.split('-')[0]+'\$-\$'+name.split('-')[-1]
            s350=Popen('grep %s table_082819.txt' %name,shell=True,stdout=PIPE).communicate()[0].split('&')[-2].split('}')[-1].strip().split('(')[0].split('n')[-1]
            p0=resn[resn['name']==nm]['p0'][0]
            tdm=(len(n)/p0)*(8.3e-6*float(dm)*df/cent_freq**3)  #bins
            ts=(len(n)/p0)*10**(-6-3.59+0.129*np.log10(float(dm))+1.02*(np.log10(float(dm)))**2-4.4*np.log10(float(cent_freq))) #bins
            err=np.sqrt(tdm**2+ts**2)/2
	    if len(nm.split('-'))>1:
		nm=nm.split('-')[0]+'$-$'+nm.split('-')[-1]
            ax[j,k].plot(n,v,lw=0.4,c='black')
            ax[j,k].set_ylim([-0.2,1.3])
            ax[j,k].set_xlim([-0.05*n[-1],1.05*n[-1]])
            ax[j,k].set_yticks([])
            ax[j,k].set_xticks([])
            ax[j,k].text(0.05*len(v),1.1*v.max(),nm,fontsize=fs_txt)
            ax[j,k].text(0.05*len(v),0.85*v.max(),"%2.1f, %s" %(dm,s350),fontsize=fs_txt)
            ax[j,k].errorbar(n[list(v).index(v.max())],-0.1,xerr=err,capsize=1,lw=0.4,c='black')
    plt.subplots_adjust(left=0.07,bottom=0.04,right=0.93,top=0.96,wspace=0,hspace=0)
    print i
    plt.savefig('plots/profiles/pulse_profiles_%s.pdf' %i)
    #plt.show()
s=0
fig,ax=plt.subplots(across,down,figsize=(6,8))
for j in range(J):
    for k in range(K):
        if s<len(profs)%(across*down):
            prof_ind=len(profs)-len(profs)%(across*down)+s
            n,v,lim=np.loadtxt(profs[prof_ind],unpack=True)
            v=v/v.max()
            nm=profs[prof_ind].split('/')[-1].split('_')[0]
            dm=resn[resn['name']==nm]['dms']
            name=nm
            if len(name.split('-'))>1:
                name=name.split('-')[0]+'\$-\$'+name.split('-')[-1]
            s350=Popen('grep %s table_082819.txt' %name,shell=True,stdout=PIPE).communicate()[0].split('&')[-2].split('}')[-1].strip().split('(')[0].split('n')[-1]
            p0=resn[resn['name']==nm]['p0']
            tdm=(len(n)/p0)*(8.3e-6*float(dm)*df/cent_freq**3)  #s 
            ts=(len(n)/p0)*10**(-6-3.59+0.129*np.log10(float(dm))+1.02*(np.log10(float(dm)))**2-4.4*np.log10(float(cent_freq))) #s
	    err=np.sqrt(tdm**2+ts**2)/2
	    if len(nm.split('-'))>1:
		nm=nm.split('-')[0]+'$-$'+nm.split('-')[-1]
            ax[j,k].plot(n,v,lw=0.4,c='black')
            ax[j,k].set_ylim([-0.2,1.3])
            ax[j,k].set_xlim([-0.05*n[-1],1.05*n[-1]])
            ax[j,k].text(0.05*len(v),1.1*v.max(),nm,fontsize=fs_txt)
            ax[j,k].text(0.05*len(v),0.85*v.max(),"%2.1f, %s" %(dm,s350),fontsize=fs_txt)
            ax[j,k].errorbar(n[list(v).index(v.max())],-0.1,xerr=err,capsize=1,lw=0.4,c='black')
        else:
            ax[j,k].set_visible(False)
        ax[j,k].set_yticks([])
        ax[j,k].set_xticks([])    
        s+=1
plt.subplots_adjust(left=0.07,bottom=0.04,right=0.93,top=0.96,wspace=0,hspace=0)
plt.savefig('plots/profiles/pulse_profiles_%s.pdf' %(i+1))
#plt.show()
