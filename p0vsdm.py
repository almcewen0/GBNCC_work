import numpy as np
import matplotlib.pyplot as plt
import subprocess as sproc
import datetime as dt
import matplotlib as mpl
import pylab
from astropy.table import Table
from glob import glob
from functions_python import wsmr,dot,s350

fontsize=14
fs_lab=14
pylab.rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
pylab.rc('xtick',labelsize=fontsize)
pylab.rc('ytick',labelsize=fontsize)
pylab.rc('text',usetex=True)

snrmin=3.8 #12.
tsys=34.4*((0.35/.408)**-2.6)  #K
g=2.
npol=2.
tint=120.
bw=80
dc=0.08
pfc=snrmin*tsys*1.13/(g*np.sqrt(npol*tint*bw))

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
full=Table.read('psr_lists/full_list.txt',format='ascii')
cond=[p['name'] not in resn['name'] and bool(p['in_surv']) for p in full]
for i in range(len(cond)):
    try:
        int(full[i]['name'][-1])
    except:
        cond[i]=False
t_nondet=full[cond]

unexp_det=[]
for p in resn:
    snew,ds=s350(3.8,p['temp'],67,0.06,1,p['dist'],1)
    if snew>p['s350m']:
	unexp_det.append(True)
    else:
	unexp_det.append(False)
p0=resn['p0']
dm=resn['dm']

unexpected=[]
rat_arr=[]
for p in t_nondet:
    if p['s350e']!='*' and p['s_curve']!='smeared out':
	rat=float(p['s350e'])/float(p['s_curve'])
	rat_arr.append(rat)
	if rat>1:
	    unexpected.append(True)
	else:
	    unexpected.append(False)
    else:
	rat_arr.append(0.)
	unexpected.append(False)
rat_arr=np.array(rat_arr)
pnd_abv=t_nondet[unexpected]['p0']
dmnd_abv=t_nondet[unexpected]['dm']
pnd=t_nondet[[not i for i in unexpected]]['p0']
dmnd=t_nondet[[not i for i in unexpected]]['dm']

plt.subplots()#figsize=(10,8))
plt.xscale('log')
plt.yscale('log')
plt.scatter(p0[[not i for i in unexp_det]],dm[[not i for i in unexp_det]],marker='+',s=8,label='Expected Detections',c='blue',alpha=0.25)
plt.scatter(p0[unexp_det],dm[unexp_det],marker='x',s=8,label='Unexpected Detections',c='blue')
plt.scatter(pnd,dmnd,marker='^',s=8,label='Expected Non-detections',color='red',alpha=0.25)
plt.scatter([],[],s=5,c='red',label='$S_{\\rm exp}$/$S_{\\rm lim}$ = 5')
plt.scatter([],[],s=20,c='red',label='$S_{\\rm exp}$/$S_{\\rm lim}$ = 20')
plt.scatter([],[],s=40,c='red',label='$S_{\\rm exp}$/$S_{\\rm lim}$ = 40')
plt.legend(loc='upper left')#,fontsize=fs_leg)
plt.scatter(pnd_abv,dmnd_abv,marker='o',c='red',s=rat_arr[unexpected])#,label='Unexpected Non-detections')

#plt.title("P0 vs. DM for Detections and Non-detections")#fontsize=fs_tit)
plt.ylabel('Dispersion Measure [pc cm$^{-3}$]',fontsize=fs_lab)
plt.xlabel('Spin Period [s]', fontsize=fs_lab)
plt.tick_params()#labelsize=fs_ticks)
dat=str((dt.datetime.now().date())).split('-')
curdat=dat[1]+dat[2]+dat[0][2:]
plt.subplots_adjust(left=0.11,right=0.98,top=0.98,bottom=0.12,wspace=0.0,hspace=0.0)
#plt.savefig('plots/p0vsDM_%s.pdf' %curdat)
plt.show()

