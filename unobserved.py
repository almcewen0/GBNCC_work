import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord

def inv(l):
    return [not i for i in l]

resp=raw_input('hist, map, or both (h/m/b)?')
p=Table.read('GBNCC_pointings.fits')
t=Table.read('GBNCC_mask_fraction.fits')
obs=Table.read('GBNCC_observed_052019.txt',format='ascii')
cond=p['obs']==0

if resp=='h' or resp=='b':
	plt.figure(figsize=(10,7))
	rfi_ra=SkyCoord(t[t['true_bw']<50]['gl'],t[t['true_bw']<50]['gb'],unit='deg',frame='galactic').transform_to('icrs').ra.value
	n,bins,patches=plt.hist([p[cond]['RAdeg'],rfi_ra],bins=150,align='mid',rwidth=0.9,stacked=True,label=['Unobserved','RFI - redo'])
	plt.xlabel('RA [deg]')
	plt.ylabel('Count')
	plt.grid(axis='y')
	plt.legend(loc='upper left')
	plt.savefig('unobserved_hist.pdf')
	plt.title('GBNCC Beams to Observe')
	plt.show()

if resp=='m' or resp=='b':
	f = plt.figure(figsize=(10,8))
	ax = f.add_subplot(111,projection='mollweide')
	pos=SkyCoord(p['RAdeg'],p['Decdeg'],unit='deg').transform_to('galactic')
	l,b=pos.l.value,pos.b.value
	inds=np.where(l>180.0)
	l[inds]-=360.0

	rfi_pos=SkyCoord(t[t['true_bw']<50]['gl'],t[t['true_bw']<50]['gb'],unit='deg',frame='galactic')
	rl,rb=rfi_pos.l.value,rfi_pos.b.value
	inds=np.where(rl>180.0)
	rl[inds]-=360.0

	deg2rad=np.pi/180
	plt.scatter(l[cond]*deg2rad,b[cond]*deg2rad,marker='.',c='blue',alpha=0.1,lw=0)
	plt.scatter(360,0,c='blue',label='Unobserved')
	plt.scatter(l[inv(cond)]*deg2rad,b[inv(cond)]*deg2rad,marker='.',c='grey',alpha=0.05,lw=0)
	plt.scatter(360,0,c='grey',label='Observed')
	plt.scatter(rl*deg2rad,rb*deg2rad,marker='.',c='red',alpha=0.5,lw=0)
	plt.scatter(360,0,c='red',label='RFI - redo')
	plt.grid()
	plt.subplots_adjust(left=0.08,bottom=0,right=0.98,top=1)
	plt.xlabel('Galactic Longitude [deg]',fontsize=15)
	plt.ylabel('Galactic Latitude [deg]',fontsize=15)
	plt.title('GBNCC Beams to Observe')
	plt.legend(loc='upper left')
	plt.savefig('unobserved_skymap.pdf')
	plt.show()


