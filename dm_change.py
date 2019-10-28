import numpy as np
from astropy.table import Table
from astropy.coordinates import SkyCoord
from functions_python import tex_error,catalog  
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


tmp=open('large_dm_change.tex','w')
j=0
t=resn[abs(resn['dms']/resn['dm']-1)>=0.05]
dmwerr=tex_error(t['dms'],t['dmserr'])     
for i,p in enumerate(t):
    nm=p['name']
    dmo=catalog(nm,'-c "JNAME DM"')
    if p['surv']!='catalog':
        if p['surv']!='gbncc':    
            nm='%s$^{\\rm %s}$' %(nm,p['surv'][0].capitalize())
        else:         
            nm=nm+'$^{\\rm B}$'           
    flg=False
    if p['surv']!='catalog' or float(dmo[4])==0:
        #print p['dm'],0.1*p['dm']
        dm_cat=str(p['dm'])
        if '.' not in dm_cat:
            dm_cat=dm_cat+'.'
        ddm_cat=0.1*p['dm']
        flg=True
    else:        
        if len(str(dmo[3]).split('.'))>1:
            #print dmo[3],float(dmo[4])*10**(-len(str(dmo[3]).split('.')[-1]))     
            ddm_cat=float(dmo[4])*10**(-len(str(dmo[3]).split('.')[-1]))
        else:    
            #print dmo[3],dmo[4]    
            ddm_cat=dmo[4]
        dm_cat=dmo[3]
    if len(nm.split('-'))>1:
        nm=nm.split('-')[0]+'$-$'+nm.split('-')[1]
    cat_werr=tex_error([dm_cat],[ddm_cat])[0]
    pad=''
    for j in range(3-len(cat_werr.split('.')[0])):
        pad=pad+'0'
    if len(pad)>0:
        cat_werr='\\phantom{%s}%s' %(pad,cat_werr)
    if flg:
        cat_werr=cat_werr+'$^{*}$'
    tmp.write('%s & %s & %s & %.3f \\\\ \n' %(nm,cat_werr,dmwerr[i],p['p0']))
tmp.close()



'''
dmwerr=tex_error(resn[cnd]['dms'],resn[cnd]['dmserr'])
for i,p in enumerate(resn[cnd]):
    nm=p['name']
    dmo=catalog(nm,'-c "JNAME DM"')
    if p['surv']!='catalog':
        if p['surv']!='gbncc':    
            nm='%s$^{\\rm %s}$' %(nm,p['surv'][0].capitalize())
        else:         
            nm=nm+'$^{\\rm B}$'           
    flg=False
    if p['surv']!='catalog' or float(dmo[4])==0:
        #print p['dm'],0.1*p['dm']
        dm_cat=str(p['dm'])
        if '.' not in dm_cat:
            dm_cat=dm_cat+'.'
        ddm_cat=0.1*p['dm']
        flg=True
    else:        
        if len(str(dmo[3]).split('.'))>1:
            #print dmo[3],float(dmo[4])*10**(-len(str(dmo[3]).split('.')[-1]))     
            ddm_cat=float(dmo[4])*10**(-len(str(dmo[3]).split('.')[-1]))
        else:    
            #print dmo[3],dmo[4]    
            ddm_cat=dmo[4]
        dm_cat=dmo[3]
    if len(nm.split('-'))>1:
        nm=nm.split('-')[0]+'$-$'+nm.split('-')[1]
    cat_werr=tex_error([dm_cat],[ddm_cat])[0]
    pad=''
    for j in range(3-len(cat_werr.split('.')[0])):
        pad=pad+'0'
    if len(pad)>0:
        cat_werr='\\phantom{%s}%s' %(pad,cat_werr)
    if flg:
        cat_werr=cat_werr+'$^{*}$'
    print "%s & %s & %s & %.3f \\\\" %(nm,cat_werr,dmwerr[i],p['p0'])
'''
