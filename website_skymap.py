import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.coordinates import SkyCoord
from matplotlib.colors import LinearSegmentedColormap    
import datetime as dt
dat=str((dt.datetime.now().date())).split('-')
curdat=dat[1]+dat[2]+dat[0][2:]

obs=Table.read('obs_forwebsite.txt',format='ascii')
pos=SkyCoord(obs['l'],obs['b'],frame='galactic',unit='deg').transform_to('icrs')
psr=Table.read('psr_lists/timing_GBNCC_111918.txt',format='ascii',names=('num','name','ra','dec','p0','p1','dm'))
psrs=SkyCoord(psr[psr['p0']>0.03]['ra'],psr[psr['p0']>0.03]['dec'],unit=('hourangle','deg'))
msp=SkyCoord(psr[psr['p0']<0.03]['ra'],psr[psr['p0']<0.03]['dec'],unit=('hourangle','deg'))  
ra,dec=np.loadtxt('psr_lists/rrats.txt',usecols=[1,2],unpack=True,dtype=str) 
rrat=SkyCoord(ra,dec,unit=('hourangle','deg')) 
deg2rad=np.pi/180.

marker=[]       
for m in obs['mjd']:
    marker.append(1+int(np.floor((m-55196)/365.)))
marker=np.array(marker)+2009     

f = plt.figure(figsize=(12,10))
ax = f.add_subplot(111,projection='mollweide')
plt.grid(b=True,linestyle='--')
ra=pos.transform_to('galactic').l.value
dec=pos.transform_to('galactic').b.value
inds=np.where(ra>180.0)
ra[inds]-=360.0
cm=LinearSegmentedColormap.from_list('Year',[(0.05,i,0.05) for i in np.linspace(0,1,11)],N=11)
rratra=rrat.transform_to('galactic').l.value
inds1=np.where(rratra>180.0)
rratra[inds1]-=360
psrra=psrs.transform_to('galactic').l.value
inds2=np.where(psrra>180.0)
psrra[inds2]-=360.0
mspra=msp.transform_to('galactic').l.value
inds3=np.where(mspra>180.0)
mspra[inds3]-=360.0
sc=plt.scatter(ra*deg2rad,dec*deg2rad,marker='.',cmap=cm,c=marker.astype(str),lw=0.2)
plt.scatter(psrra*deg2rad,psrs.transform_to('galactic').b.value*deg2rad,edgecolor='black',c='white',marker='o',linewidth=1,s=30,label='Normal Pulsars')
plt.scatter(rratra*deg2rad,rrat.transform_to('galactic').b.value*deg2rad,marker='^',edgecolor='black',c='white',label='RRATs',s=90)
plt.scatter(mspra*deg2rad,msp.transform_to('galactic').b.value*deg2rad,edgecolor='black',c='white',marker='o',linewidth=1,s=90,label='MSPs')
plt.legend(loc='lower right',fontsize=12)
plt.colorbar(sc,shrink=0.5,ticks=[2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019],spacing='uniform')
plt.title('GBNCC Progress: %i beams observed' %len(ra),fontsize=15)
plt.ylabel('Galactic Latitude [deg]',fontsize=15)
plt.xlabel('Galactic Longitude [deg]',fontsize=15) 
plt.savefig('plots/gbncc_skymap_glgb_%s.png' %curdat)
#plt.show()



f = plt.figure(figsize=(12,10))
ax = f.add_subplot(111,projection='mollweide')
plt.grid(b=True,linestyle='--')
ra=pos.ra.value
dec=pos.dec.value
inds=np.where(ra>180.0)
ra[inds]-=360.0
cm=LinearSegmentedColormap.from_list('Year',[(0.05,i,0.05) for i in np.linspace(0,1,11)],N=11)
rratra=rrat.ra.value
inds1=np.where(rratra>180.0)
rratra[inds1]-=360
psrra=psrs.ra.value
inds2=np.where(psrra>180.0)
psrra[inds2]-=360.0
mspra=msp.ra.value
inds3=np.where(mspra>180.0)
mspra[inds3]-=360.0
sc=plt.scatter(ra*deg2rad,dec*deg2rad,marker='.',cmap=cm,c=marker.astype(str),lw=0.2)
plt.scatter(psrra*deg2rad,psrs.dec.value*deg2rad,edgecolor='black',c='white',marker='o',linewidth=1,s=30,label='Normal Pulsars')
plt.scatter(rratra*deg2rad,rrat.dec.value*deg2rad,marker='^',edgecolor='black',c='white',label='RRATs',s=90)
plt.scatter(mspra*deg2rad,msp.dec.value*deg2rad,edgecolor='black',c='white',marker='o',linewidth=1,s=90,label='MSPs')
plt.legend(loc='lower right',fontsize=12)
plt.colorbar(sc,shrink=0.5,ticks=[2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019],spacing='uniform')
plt.title('GBNCC Progress: %i beams observed' %len(ra),fontsize=15)
plt.ylabel('Declination [deg]',fontsize=15)
plt.xlabel('Right Ascension [deg]',fontsize=15) 
plt.savefig('plots/gbncc_skymap_radec_%s.png' %curdat)
#plt.show()

