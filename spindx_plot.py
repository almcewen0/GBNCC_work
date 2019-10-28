import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table,Column
from functions_python import dot
import scipy.optimize as opt
from glob import glob


def line(x,m,b):
    return m*x+b


sdata=Table.read('catalog_fluxes_psrsinsurv.txt',format='ascii')
nms=list(sdata[0]) 
for i in range(len(sdata.colnames)):
    sdata[sdata.colnames[i]].name=nms[i]
sdata.remove_row(0)
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
tmp=open('spindx_measurements.txt','w')

spindx_meas=[] 
dspindx_m=[] 
for p in sdata: 
    if p['name'] in np.loadtxt('broken_spindx_psrs.txt',dtype=str):
        spindx_meas.append('*')
        dspindx_m.append('*')
        continue
    x=[] 
    y=[] 
    yerr=[] 
    spindx=p['spindx'] 
    for i in range(5,39,2): 
        if float(nms[i])<350 and float(nms[i+2])>350: 
            x.append(350.) 
            y.append(resn[resn['name']==p['name']]['s350m'][0]) 
            yerr.append(resn[resn['name']==p['name']]['ds'][0]) 
        if p[i]!='*': 
            x.append(float(nms[i])) 
            y.append(float(p[i])) 
            yerr.append(float(p[i+1])) 
    x=np.array(x) 
    y=np.array(y) 
    #delta=np.sum(1./y.std()**2)*np.sum(x**2/y.std()**2)-np.sum(x/y.std()**2)**2 
    if len(x)<3:
	tmp.write('%s * *\n' %p['name'])
        spindx_meas.append('*')
        dspindx_m.append('*')
        continue
        #f,cov=opt.curve_fit(line,np.log10(x),np.log10(y),sigma=yerr/y,p0=[-2,10]) 
    else:
        f,cov=opt.curve_fit(line,np.log10(x),np.log10(y),p0=[-2,10]) 
    f=np.poly1d(f)
    f_err=np.sqrt(np.diag(cov))  #np.sqrt(sum((np.log10(x)-np.log10(x).mean())**2)/((len(x)-2)*sum((np.log10(x)-f(np.log10(x)))**2))) #(1./delta)*(np.sum(1./y.std()**2)*np.sum(x*y/y.std()**2)-np.sum(x/y.std()**2)*np.sum(y/y.std()**2)) 
    spindx_meas.append(f[1])
    if f[1]!=p['spindx_meas']:
        p['spindx_meas']=float(f[1])
        p['dspindx_m']=float(f_err[0])
    if str(f_err[0])=='inf':
        f_err=[((yerr[-1]/y[-1])+(yerr[0]/y[0]))/(np.log10(x[-1])-np.log10(x[0]))]
    dspindx_m.append(f_err[0]) 
    plt.plot(np.log10(x),np.log10(y),'+',c='red') 
    if spindx!='*':
        x0=(np.log10(x[-1])+np.log10(x[0]))/2 
        b2=f(x0)-float(spindx)*x0 
        #plt.plot(np.log10(x),b2+float(spindx)*(np.log10(x)),lw=1,label=r'$\alpha_{pub} = %2.2f \pm %2.2f$' %(float(spindx),float(p['dspindx'])),ls='--',c='orange') 
    tmp.write('%s %f %f\n' %(p['name'],spindx_meas[-1],dspindx_m[-1]))
    plt.plot(np.log10(x),f(np.log10(x)),lw=1,label=r'$\alpha_{BF} = %2.2f \pm %2.2f$' %(f[1],f_err[0]),ls='-.',c='blue') 
    yerr=np.array(yerr)/np.array(y) 
    plt.errorbar(np.log10(x),np.log10(y),yerr=yerr,lw=0,elinewidth=1,capsize=2,label='Flux Measurements',c='red') 
    plt.axvline(x=np.log10(350),ls='--',lw=0.5,c='orange',label='350 MHz') 
    plt.ylabel('Measured Flux') 
    plt.xlabel('Observing Frequency') 
    plt.xticks([np.log10(x[0]),np.log10(x[0])+(np.log10(x[-1])-np.log10(x[0]))/3,np.log10(x[0])+2*(np.log10(x[-1])-np.log10(x[0]))/3,np.log10(x[-1])],['%s' %int(x[0]),'%s'%int(10**(np.log10(x[0])+(np.log10(x[-1])-np.log10(x[0]))/3)),'%s'%int(10**(np.log10(x[0])+2*(np.log10(x[-1])-np.log10(x[0]))/3)),'%s' %int(x[-1])]) 
    plt.yticks(plt.yticks()[0][0::2],['10$^{%1.0f}$' %v for v in plt.yticks()[0][0::2].astype(float)]) 
    plt.legend() 
    plt.suptitle(p['name']) 
    plt.savefig('/Users/aemcewen/GBNCC_paper/plots/spindx/%s.png' %p['name']) 
    #plt.show() 
    plt.clf() 
    #plt.close()
tmp.close()

#sdata.remove_column('spindx_meas')
#sdata.remove_column('dspindx_m')
#spindx_meas=Column(spindx_meas,name='spindx_meas')
#dspindx_m=Column(dspindx_m,name='dspindx_m')
#sdata.add_column(spindx_meas)
#sdata.add_column(dspindx_m)
#sdata.write('catalog_fluxes_psrsinsurv.txt',format='ascii',overwrite=True)    
