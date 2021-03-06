import numpy as np
import datetime as dt
from subprocess import call,Popen,PIPE
from glob import glob
from astropy.table import Table

resn=Table.read('../master_results_081519.txt',format='ascii')
bms=Table.read('../GBNCC_pointings.fits')
today=dt.date.today()
date_str=str(today.year)+'-'+str(today.month)+'-'+str(today.day)

def inv(l):
    return [not i for i in l]


for p in resn:
    nf,vf,limf=np.loadtxt(glob('../presto_plots/master_112818/good_profs_withwidth/*%s*width' %p['name'])[0],unpack=True)
    vf=vf-vf[limf.astype(bool)].mean()
    bmnm=str(p['bm'])
    while len(bmnm)<5:
	bmnm='0'+bmnm
    bm=bms[bms['pointing']=='GBNCC'+bmnm]
    sn=Popen("grep %s ../table_091819.txt | awk -F'&' '{print $7}'" %p['name'],shell=True,stdout=PIPE).communicate()[0].split('\n')[0].split()[0].strip('\\phn')
    s=Popen("grep %s ../table_091819.txt | awk -F'&' '{print $8}'" %p['name'],shell=True,stdout=PIPE).communicate()[0].split('\n')[0].split()[0].strip('\\phn')
    s_val=float(s.split('(')[0])
    gamma=s_val/sum(vf)
    tmp=open('%s_gbnccprof.txt' %p['name'],'w')
    tmp.write('# PSR %s \n' %p['name'])
    tmp.write('# P0 = %s s \n' %p['p0'])
    if p['surv'] != 'catalog':
        if p['surv']=='lofar':
            tmp.write('# PSR from LOTAAS survey discovery data \n')
        else:
            tmp.write('# PSR from %s survey discovery data \n' %p['surv'].upper())
    tmp.write('# Brightest detection in GBNCC Beam %s on MJD %s \n' %(p['bm'],p['mjd']))
    tmp.write('# Beam Position (RA/Dec): %s %s \n' %(bm['RA'][0],bm['Dec'][0]))
    tmp.write('# Position offset from pulsar = %s degrees \n' %p['dist'])
    tmp.write('# Tsys = %s K \n' %p['temp'])
    tmp.write('# Bandwidth after flagging and removing rolloff = %s MHz \n' %p['bw'])
    tmp.write('# Observation length = 120 s \n')
    tmp.write('# Number of on-pulse bins = %i \n' %int(p['on-pulse']))
    tmp.write('# GBT gain at 350 MHz = 2 K/Jy \n')
    tmp.write('# Beta = 1.13 \n')
    tmp.write('# S/N measured from profile = %s \n' %sn)
    tmp.write('# Flux Density at 350MHz = %s mJy \n' %s)
    tmp.write('# Profile scaled such that S350 = sum(Pn) where Pn are the profile bins using factor gamma = %s mJy/bin \n' %gamma)
    tmp.write('# Profile generated by Alex McEwen on %s using prof_maker.py and functions in functions_python on github.com/almcewen0/GBNCC_work \n\n' %date_str)
    tmp.write('# Bin Intensity On-pulse \n')
    for n,v,l in zip(nf,vf,limf):
        tmp.write('%i %s %s \n' %(int(n),v*gamma,int(not l)))
    tmp.close()
