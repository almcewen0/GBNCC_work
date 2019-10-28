import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
import scipy.optimize as opt
from functions_python import tex_error
from glob import glob 

def line(x,m,b):
    return m*x+b

#figsize=(10,8)
#import matplotlib as mpl
#fs_lab=20
#fs_leg=14
#fs_tit=20
#fs_ticks=12
#mpl.rc('font',family='Arial')

import pylab

fontsize=14
fs_lab=14
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
sdata=Table.read('catalog_fluxes_psrsinsurv.txt',format='ascii')
nms=list(sdata[0]) 
for i in range(len(sdata.colnames)):
    sdata[sdata.colnames[i]].name=nms[i]
sdata.remove_row(0)

fig,ax=plt.subplots(2,2)#,figsize=figsize)
j,k=0,0
name,a_l,da_l,a_h,da_h,brk_arr=[],[],[],[],[],[]
print "PSR & \\alpha_{\\rm{l}} & \\alpha_{\\rm{h}} & Break Frequency [MHz]\\\ "
for nm in np.loadtxt('broken_spindx_psrs.txt',dtype=str):
    p=sdata[sdata['name']==nm][0]
    x=[]
    y=[]
    yerr=[]
    spindx=p['spindx']
    for i in range(5,40,2):
        if float(nms[i])<350 and float(nms[i+2])>350:
            x.append(350.)
            y.append(resn[resn['name']==p['name']]['s350m'][0])
            yerr.append(resn[resn['name']==p['name']]['ds'][0])
        if p[i]!='*':
            x.append(float(nms[i]))
            y.append(float(p[i]))
            yerr.append(float(p[i+1]))
    if len(x)<3:
        continue
    x=np.array(x)
    y=np.array(y)
    yerr=np.array(yerr/y)
    x0=np.log10(x)
    y0=np.log10(y)
    std_arr=[]
    n=[]
    f=np.poly1d(np.polyfit(x0,y0,1))
    stde_full=(x0-f(x0)).std()/np.sqrt(len(x0)-1)
    dstd=[0.]
    for N in range(2,len(x)):
        x0=np.log10(x[:N])
        y0=np.log10(y[:N])
        tmp=np.poly1d(np.polyfit(x0,y0,1))
        std_arr.append((y0-tmp(x0)).std())
        n.append(N)
        if N>2:
            dstd.append(std_arr[-1]-std_arr[-2])
    N=n[list(dstd==np.array(dstd).max()).index(1)]-1
    xf=np.log10(x)
    x0=np.log10(x[:N])
    y0=np.log10(y[:N])
    x1=np.log10(x[N-1:]) 
    y1=np.log10(y[N-1:]) 
    if len(x0)>2:
        tmp0,cov0=opt.curve_fit(line,x0,y0,p0=[1,10])
        err0=np.sqrt(np.diag(cov0))[0]
    else:
        m=(y0[1]-y0[0])/(x0[1]-x0[0])
        b=y0[0]-m*x0[0]
        tmp0=[m,b]
        smin=((y0[1]+yerr[1])-(y0[0]-yerr[0]))/(x0[1]-x0[0])
        smax=((y0[1]-yerr[1])-(y0[0]+yerr[0]))/(x0[1]-x0[0])
        err0=abs(smax-smin)/2.
    if len(x1)>2:
        tmp1,cov1=opt.curve_fit(line,x1,y1,p0=[-1,20])
        err1=np.sqrt(np.diag(cov1))[0]
    else:
        m=(y1[1]-y1[0])/(x1[1]-x1[0])
        b=y1[0]-m*x1[0]
        tmp1=[m,b]
        smin=((y1[1]+yerr[-1])-(y1[0]-yerr[-2]))/(x1[1]-x1[0])
        smax=((y1[1]-yerr[-2])-(y1[0]+yerr[-2]))/(x1[1]-x1[0])
        err1=abs(smax-smin)/2.
    tmp0=np.poly1d(tmp0)
    tmp1=np.poly1d(tmp1)
    xt=np.arange(x[0],x[-1],0.5)
    for tp in np.log10(xt):
        if abs(tmp0(tp)-tmp1(tp))<1e-2:
            brk=tp
            break
    xl=list(x0)
    if xl[0]>brk:
        xl.insert(0,brk)
    else:
        for i in range(1,len(xl)):
            if xl[i]>brk and xl[i-1]<brk:
                xl.insert(i,brk)
    xl2=list(x1)
    if xl2[0]>brk:
        xl2.insert(0,brk)
    else:
        for i in range(1,len(xl2)):
            if xl2[i]>brk and xl2[i-1]<brk:
                xl2.insert(i,brk)
    if brk>xl[-1]:
	xl.append(brk) 
    brk_arr.append(brk)
    ax[k,j].errorbar(np.log10(x),np.log10(y),yerr=yerr,lw=0,elinewidth=1,capsize=2,c='blue',label='Measurements')
    #ax[k,j].plot(xf,f(xf),ls='--',linewidth=1,label=r'$\alpha_{\rm{BF}}$')
    ax[k,j].plot(np.array(xl)[xl<=brk],tmp0(np.array(xl)[xl<=brk]),linewidth=1,c='orange',label=r'$\alpha_{\rm{l}}$')
    ax[k,j].plot(np.array(xl2)[xl2>=brk],tmp1(np.array(xl2)[xl2>=brk]),linewidth=1,c='orange',label=r'$\alpha_{\rm{h}}$')
    ax[k,j].axvline(brk,linestyle='--',linewidth=1,c='red',label='Break Frequency')
    ax[k,j].set_ylim([-1,3])
    ax[k,j].set_xticks([np.log10(x[0]),np.log10(x[0])+(np.log10(x[-1])-np.log10(x[0]))/2,np.log10(x[-1])])
    ax[k,j].set_xticklabels(['%s' %int(x[0]),'%s'%int(10**(np.log10(x[0])+(np.log10(x[-1])-np.log10(x[0]))/2)),'%s' %int(x[-1])])#,fontsize=fs_ticks)
    ax[k,j].set_yticks([-1,0,1,2,3])
    ax[k,j].set_yticklabels(['10$^{%1.0f}$' %v for v in ax[k,j].get_yticks().astype(float)])#,fontsize=fs_ticks)
    ax[k,j].text(np.linspace(xf.min(),xf.max(),10)[6],2,nm)#,fontsize=fs_leg)
    #ax[k,j].set_title(p['name'])#,fontsize=fs_tit)
    if j==0:
        ax[k,j].set_ylabel('Flux Density [mJy]',fontsize=fs_lab)
    if k==1:
        ax[k,j].set_xlabel('Frequency [MHz]',fontsize=fs_lab)
    ax[k,j].legend(loc='lower left')#fontsize=fs_leg,
    j+=1
    if j%2==0:
        k+=1
        j-=2
        if k%2==0:
            k-=1
    name.append(nm)
    a_l.append(tmp0[1])
    da_l.append(err0)
    a_h.append(tmp1[1])
    da_h.append(err1)
a_l=tex_error(a_l,da_l)
a_h=tex_error(a_h,da_h)
for i in range(len(name)):
    print "%s & %s & $-$%s & %2.0f \\\ " %(name[i],a_l[i],a_h[i],10**brk_arr[i]) 
plt.tight_layout()
plt.savefig('plots/broken_spindx.pdf')   
#plt.show()

