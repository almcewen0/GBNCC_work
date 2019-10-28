import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import datetime as dt
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units as u
dat=str((dt.datetime.now().date())).split('-')
curdat=dat[1]+dat[2]+dat[0][2:]
from glob import glob

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
pos=SkyCoord(resn['rajd'],resn['decjd'],unit='deg').transform_to('galactic')
gb,gl=np.array(pos.b/u.deg),np.array(pos.l/u.deg)
s=resn['surv']
inds = np.where(gl > 180.0)
gl[inds] -= 360.0
deg2rad = np.pi/180.

full=Table.read('psr_lists/full_list.txt',format='ascii')
nondet=[l['name'] not in resn['name'] and l['in_surv']==1 for l in full]
posnd=SkyCoord(full[nondet]['gl'],full[nondet]['gb'],unit='deg',frame='galactic')
glnd=posnd.l.value
gbnd=posnd.b.value
inds=np.where(glnd>180.0)
glnd[inds]-=360.

f = plt.figure()#figsize=figsize)
ax = f.add_subplot(111,projection='mollweide')
outfile = 'plots/psr_map_Gal_%s.pdf' %curdat
gbnccl,gbnccb=np.loadtxt('GBNCC_observed_glgb.dat',usecols=[1,2],unpack=True,skiprows=1)
inds = np.where(gbnccl > 180.0)
gbnccl[inds] -= 360.0


plt.scatter(gbnccl*deg2rad,gbnccb*deg2rad,marker='.',c='grey',lw=0.2,alpha=0.025)
plt.scatter(glnd*deg2rad,gbnd*deg2rad,edgecolor='',marker='x',s=30,c='blue',label='Nondetections',lw=1.5)
plt.scatter(gl[s!='catalog']*deg2rad,gb[s!='catalog']*deg2rad,edgecolor='',marker='^',s=30,c='red',label='Discovery Data',lw=1.5)
plt.scatter(gl[s=='catalog']*deg2rad,gb[s=='catalog']*deg2rad,edgecolor='',marker='+',s=30,c='green',label='ATNF catalog',lw=1.5) 
print "%i from catalog" %sum(s=='catalog')
print "%i not from catalog" %sum(s!='catalog')

plt.grid()
#plt.title('Pulsars detected in GBNCC Survey')#,fontsize=fs_tit)
plt.legend(loc='lower right')#,fontsize=fs_leg)
plt.xlabel('Galactic Longitude [deg]',fontsize=fs_lab)
plt.ylabel('Galactic Latitude [deg]',fontsize=fs_lab)
plt.tick_params()#labelsize=fs_ticks)
plt.savefig(outfile)	
#plt.show()
