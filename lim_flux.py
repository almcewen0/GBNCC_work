import numpy as np 
from astropy.table import Table
import matplotlib.pyplot as plt
from glob import glob

#figsize=(10,8)
#import matplotlib as mpl
#fs_leg=16
#fs_tit=22
#fs_ticks=14
#mpl.rc('font',family='Arial')

import pylab

fontsize=14
fs_lab=20
pylab.rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
pylab.rc('xtick',labelsize=fontsize)
pylab.rc('ytick',labelsize=fontsize)
pylab.rc('text',usetex=True)



resp=raw_input('hist, map, or both (h/m/b)?')
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
t = Table.read('GBNCC_mask_fraction.fits')


if resp=='h' or resp=='b':
	SN = 3.8 	# limiting S/N
	G = 2. 		# K/Jy
	Npol = 2. 	# 2 polarizations
	Tobs = 120. 	# 120-sec scans
	beta = 1.3 	# instrumental correction factor
	f = 1 		# beam factor

	duty = np.average(resn['prof_len']/resn['on-pulse'])-1 							# p/w-1
	flux = t['flux_lim'] #SN * np.array(tsys) * beta / (G * np.sqrt(duty * Npol * Tobs * np.array(frac)) )  # radiometer equation
	cond=flux.astype(str)!='inf' 										# some beams were completely masked, so flux_lim=infinity
	print("Minimum S/N: %i" %SN)
	print("Minimum limiting flux density: %f mJy" % min(flux[cond]))
	print("Maximum limiting flux density: %f mJy" % max(flux[cond]))
	print("Mean limiting flux density: %f mJy" % np.average(flux[cond]))
	print("Median limiting flux density: %f mJy" % np.median(flux[cond]))

	plt.figure()#figsize=figsize)
	plt.xlim([0.35,13])
	plt.hist(flux[cond], bins=1000, color = 'blue', cumulative = True, normed=True, histtype = 'step')
	plt.ylabel('Fractional Number of GBNCC Pointings',fontsize=fs_lab)
	plt.xlabel('Limiting Flux Density [mJy]',fontsize=fs_lab)
	#plt.title('Limiting Flux Density using S/N of 10 for all GBNCC pointings')
	plt.tick_params()#labelsize=fs_ticks)
	plt.xscale('log')

	txtloc=np.linspace(0.7,0.8,4)
	#plt.text(5,txtloc[3],"Minimum:      %2.2f mJy" % min(flux[cond]))#,fontsize=fs_leg)
	#plt.text(5,txtloc[2],"Maximum:   %2.2f mJy" % max(flux[cond]))#,fontsize=fs_leg)
	#plt.text(5,txtloc[1],"Mean:            %2.2f mJy" % np.average(flux[cond]))#,fontsize=fs_leg)
	#plt.text(5,txtloc[0],"Median:         %2.2f mJy" % np.median(flux[cond]))#,fontsize=fs_leg)
	plt.subplots_adjust(left=0.11,bottom=0.12,right=0.98,top=0.98)
	plt.savefig('plots/GBNCC_fluxdensity_limit_hist.pdf')
	#plt.show()

if resp=='m' or resp=='b':
	f = plt.figure()#figsize=figsize)
	ax = f.add_subplot(111,projection='mollweide')
	l=t['gl']
	b=t['gb']
	inds=np.where(l>180.0)
	l[inds]-=360.0
	deg2rad=np.pi/180
	cm=plt.cm.get_cmap('plasma')
	plt.scatter(l*deg2rad,b*deg2rad,marker='.',cmap=cm,c=np.log10(t['flux_lim']),alpha=0.1,lw=0.001,vmax=0.5)
	sc=plt.scatter((200,200),(0,0),cmap=cm,c=[t['flux_lim'].max(),t['flux_lim'].min()],vmax=10**0.5)
	cb=plt.colorbar(sc,shrink=0.5)
	cb.set_label('Limiting Flux Density [mJy]',fontsize=fs_lab-6)
	cb.ax.tick_params(labelsize=fs_lab-6)
	plt.grid()
	plt.subplots_adjust(left=0.08,bottom=0,right=1,top=1)
	plt.tick_params()#labelsize=fs_ticks)
	plt.xlabel('Galactic Longitude [deg]',fontsize=fs_lab-6)
	plt.ylabel('Galactic Latitude [deg]',fontsize=fs_lab-6)
	plt.subplots_adjust(left=0.14,bottom=0,right=1,top=1)
	plt.savefig('plots/limiting_flux_skymap.pdf')
	#plt.show()


