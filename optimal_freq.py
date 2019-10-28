import numpy as np 
import frequencyoptimizer as fo
from astropy.table import Table
from glob import glob
from functions_python import dot

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

galnoise = fo.GalacticNoise()
telnoise = fo.TelescopeNoise(gain=2.0,T_const=30)

print ""
print "using values that were manually taken from NE2001 website - update if S/Ns change!"
print ""

d=   [1.737,   1.520,   0.503,    0.747,    1.469,   2.935,    2.597,   3.307,   0.913,   1.053]      #distance (kpc)
dnud=[0.7e-3,  1.25e-3, 21.77e-3, 7.81e-3,  6.39e-3, 1.75e-3,  0.12e-3, 1.29e-3, 3.3e-3,  5.54e-3]   #scintillation bandwidth (GHz)
dtd= [280.0,   350.6,   842.8,    615.5,    780.4,   576.7,    141.3,   526.4,   442.5,   615.1]     #scintillation timescale (s)

NCHAN = 50
NSTEPS = 50
freq_bs=[]
freq_cs=[]
mins=[]
psr_list=resn[dot(dot(resn['p0']<0.03,resn['surv']=='gbncc'),resn['snm']<=7*resn['snm'].min())]
res_tbl=Table(names=('name','p0','dm','opt_freq_GHz','opt_bw_GHz','sig_toa_min_us','sopt_mJy'),dtype=('S10',float,float,float,float,float,float))
for i,p in enumerate(psr_list):
    if p['spindx_m']!='*':
        alpha=abs(float(p['spindx_m']))
    elif p['spindx']!='*':
        alpha=abs(float(p['spindx']))
    else:
        alpha=1.6
    taud=10**(-6.46+0.154*np.log10(p['dms'])+1.07*np.log10(p['dms'])**2-3.86*np.log10(1))
    dt=724*d[i]**(-0.5)*(taud/0.0001)**(-0.5)
    psrnoise = fo.PulsarNoise(p['name'],alpha=alpha,I_0=p['s350m']*(1000./350)**(-alpha),taud=taud,DM=p['dms'],D=d[i],dnud=dnud[i],tauvar=12.2e-3,dtd=dt,Weffs=1000*p['p0']*p['on-pulse']/p['prof_len'],W50s=p['w50']*1000.,sigma_Js=np.zeros(NCHAN)+0.066,P=p['p0']*1000.)  #jitter upper limit

    freqopt = fo.FrequencyOptimizer(psrnoise,galnoise,telnoise,numin=0.1,numax=2,nchan=NCHAN,log=True,nsteps=NSTEPS)
    freqopt.calc()
    freqopt.plot("test",doshow=False,minimum='k*',save=False,points=(1.3,1.2,'ko'))
    #freqopt.save('test')
    checkdata=np.log10(freqopt.sigmas)
    flatdata=checkdata.flatten()
    inds=np.where((~np.isnan(flatdata))&~(np.isinf(flatdata)))[0]
    MIN = np.min(flatdata[inds])
    INDC,INDB = np.where(checkdata==MIN)
    INDC,INDB = INDC[0],INDB[0]
    sopt=p['s350m']*(freqopt.Bs[INDB]/0.35)**(-alpha)
    res_tbl.add_row([p['name'],p['p0'],p['dm'],freqopt.Bs[INDB],freqopt.Cs[INDC],10**MIN,sopt])
res_tbl.write('optimal_frequency.txt',format='ascii',overwrite=True)
print res_tbl
