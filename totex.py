import subprocess as sproc
import numpy as np
from astropy.table import Table
from functions_python import tex_error
from glob import glob

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
dmerr=tex_error(resn['dms'],resn['dmserr'])
snerr=tex_error(resn['snm'],resn['dsnm'])
s350err=tex_error(resn['s350m'],resn['ds'])
spindxerr=tex_error(resn['spindx_m'],resn['dspindx_m'])

#      nm    dm  dmerr  mjd   dist    w50     w10    snm dsnm       sm  dsm  spindx
s=0
null=np.loadtxt('nuller_psrnames.txt',dtype='str')
for p in resn:
    nm=p['name']
    if len(nm.split('-'))>1:
	nm=nm.split('-')[0]+'$-$'+nm.split('-')[1]
    if p['name'] in null:
	nm=nm+'*'
    if p['surv']!='catalog':
	if p['surv']=='gbncc':
	    nm=nm+'$^{\\rm{B}}$'
	elif p['surv']=='aodrift':
	    nm=nm+'$^{\\rm{A}}$'
        elif p['surv']=='gbt350':
            nm=nm+'$^{\\rm{G}}$'
        elif p['surv']=='lofar':
            nm=nm+'$^{\\rm{L}}$'
        elif p['surv']=='htru_mlt':
            nm=nm+'$^{\\rm{ML}}$'
        elif p['surv']=='htru_hlt':
            nm=nm+'$^{\\rm{HL}}$'
    nm='\object[PSR %s]{%s}' %(p['name'],nm)
    if p['w10']!='*':
	w10="%s" %int(p['w10'])
	if w10=='0':
	    w10="%1.1f" %float(p['w10'])
	if len(w10.split('.')[0])<3:
	    diff=3-len(w10.split('.')[0])
	    for i in range(diff):
		w10='\phn '+w10
    else:
	w10='*'
    w50="%s" %int(p['w50'])
    if w50=='0':
	w50="%1.1f" %float(p['w50'])
    if len(w50.split('.')[0])<3:
	diff=3-len(w50.split('.')[0])
	for i in range(diff):
	    w50='\phn '+w50
    bm=str(p['bm'])
    if len(bm)<5:
	diff=5-len(bm)
	for i in range(diff):
	    bm='0'+bm
    mjd=p['mjd']
    dist=float(p['dist'])
    numdet=p['numdet']
    if p['nosearch']==0:
        dmwerr=dmerr[s] #sproc.Popen("grep %s dms_witherr_022519.txt | awk '{print $2}'" %l[0],shell=True,stdout=sproc.PIPE).communicate()[0].split('\n')[0]
    else:
	dmwerr=str(p['dm'])+'$^{\ddagger}$'
    snwerr=snerr[s] #sproc.Popen("grep %s snm_witherr_022519.txt | awk '{print $2}'" %l[0],shell=True,stdout=sproc.PIPE).communicate()[0].split('\n')[0]
    s350werr=s350err[s] #sproc.Popen("grep %s s350_witherr_022519.txt | awk '{print $2}'" %l[0],shell=True,stdout=sproc.PIPE).communicate()[0].split('\n')[0]
    if p['name'] in np.loadtxt('broken_spindx_psrs.txt',dtype=str):
	spindx='\phn\phn\dag'
    else:
	if spindxerr[s]!='$-$':
	    if len(spindxerr[s].split('-'))>1:
    	        spindx=spindxerr[s].split('-')[0]+'$-$'+spindxerr[s].split('-')[1]
	    else:
	        spindx=spindxerr[s]
	else:
	    spindx='\\phn\\phn'+spindxerr[s]
    s+=1
    print "%s & %s & %s & %2.3f & %s & %s & %s & %s & %s \\\ " %(nm,dmwerr,mjd,dist,w50,w10,snwerr,s350werr,spindx)    
