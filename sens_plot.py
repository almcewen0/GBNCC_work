import numpy as np
import matplotlib.pyplot as plt
import pylab
import datetime as dt
from astropy.table import Table
from glob import glob
from functions_python import wsmr,dot

fontsize=14
fs_lab=20
pylab.rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
pylab.rc('xtick',labelsize=fontsize)
pylab.rc('ytick',labelsize=fontsize)
pylab.rc('text',usetex=True)

snrmin=6
g=2.
npol=2.
tint=120.
bw=67
dc=0.06

res=glob('master_results_*.txt')
if len(res)>1:
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

t=Table.read(res,format='ascii')
full=Table.read('psr_lists/full_list.txt',format='ascii')
cond=[p not in t['name'] for p in full['name']]

step=10
dms,temps=[],[]
for dm in range(10,440,step):
    cnd=dot(t['dms']<dm,t['dms']>dm-step)
    if sum(cnd)>0:
        dms.append(dm)
        temps.append(t[cnd]['temp'].mean())
temp_avg=np.poly1d(np.polyfit(dms,temps,2))

def pfc_calc(tsys,snrmin=6.,beta=1.1,g=2,npol=2,tint=120,bw=67):
    return snrmin*beta*tsys/(g*np.sqrt(npol*tint*bw))


t_nondet=full[cond]
flxminc=t['s350m']==t['s350m'].min()
pfc=pfc_calc(tsys=t[flxminc]['temp'])
offset_width=wsmr(t[flxminc]['p0']*t[flxminc]['on-pulse']/t[flxminc]['prof_len'],t[flxminc]['dms'])
offset=t[flxminc]['s350m']/(pfc*np.sqrt(offset_width/(t[flxminc]['p0']-offset_width)))
pfc=pfc*offset

x=np.arange(0.00005,1.5,0.00001)
fig,ax=plt.subplots()#figsize=figsize)
cm=plt.cm.get_cmap('RdYlBu')
for j in [20,50,100,150,200,300]:
    pfc=pfc_calc(tsys=temp_avg(j))*offset
    ax.plot(x,pfc*np.sqrt(wsmr(x*dc,j)/(x-wsmr(x*dc,j))),color=cm(j),label='%s pc cm$^{-3}$' %j)
ax.set_xscale('log',basex=10.)
ax.set_yscale('log',basey=10.)
cond2=t_nondet['s350e']!='*'
ax.scatter(t_nondet['p0'][cond2],t_nondet['s350e'][cond2].astype(float),c=t_nondet[cond2]['dm'],cmap=cm,vmax=300,lw=0.1,marker='^',label='Nondetections')
cond=t['s350m']!='*'
sc=ax.scatter(t['p0'][cond],t['s350m'][cond],c=t['dm'][cond],cmap=cm,marker='+',vmax=300,lw=1,label="Detections")
cb=plt.colorbar(sc)
cb.ax.tick_params()#labelsize=fs_ticks)
cb.set_label(label='Dispersion Measure [pc cm$^{-3}$]',fontsize=fs_lab)
ax.set_xlim(x.min())
ax.legend()#fontsize=fs_leg)
ax.tick_params()#labelsize=fs_ticks)
#plt.suptitle('Sensitivity in GBNCC Survey')#,fontsize=fs_tit)
ax.set_xlabel('Spin Period [s]',fontsize=fs_lab)
ax.set_ylabel('Flux Density [mJy]',fontsize=fs_lab)
dat=str((dt.datetime.now().date())).split('-')
curdat=dat[1]+dat[2]+dat[0][2:]
plt.subplots_adjust(left=0.13,right=0.98,top=0.98,bottom=0.13,wspace=0.0,hspace=0.0)

print "Minimum S/N = %2.3f" %(snrmin*offset)
print "Duty Cycle = %i percent" %(int(dc*100))
plt.savefig('plots/GBNCC_senscurve_%s.pdf' %curdat)
#plt.show()


#for i in range(6,15,2):
#	print i
#	sens_plot(i,0.08)

