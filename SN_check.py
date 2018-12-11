import numpy as np
from glob import glob 
import subprocess as sproc
from astropy.table import Table
import math

psr=np.loadtxt('SN_psrcheck.txt',dtype='str',usecols=[0])
beam=np.loadtxt('SN_psrcheck.txt',dtype='str',usecols=[1])
t=Table.read('psrcat_data_withtemp.fits')
centfreq=350.    #MHz
df=0.0244140625  #MHz
P=50.            #bins
width_flag=False
print "Name \t\t beam \t bw \t poln \t sep \t spindx  s350 \t temp \t\t W \t snr_data \t w50b \t snr_catW \t weffb \t snr_smearW "
out=open('SN_check_output.txt','w')
out.write("Note: ! indicates lack of catalog value for W50, so 6% of the period is used. Also, *s in the snr_smearW column indicate cases where the width is smeared out beyond the period width.\n")
out.write("\n")
out.write("Name \t\t beam \t bw \t poln \t sep \t spindx  s350 \t temp \t\t W \t snr_data \t w50b \t snr_catW \t weffb \t snr_smearW \n")
for i in range(len(psr)):
    if i>0 and psr[i]==psr[i-1]:
	chee=0
    else:
	    if t[t['name']==psr[i]]['SPINDX'][0] != '*':
		spindx=float(t[t['name']==psr[i]]['SPINDX'][0])
	    else:
		spindx=-1.7
	    if t[t['name']==psr[i]]['S400'][0] != '*':
		s350=float(t[t['name']==psr[i]]['S400'][0])*(0.35/0.4)**spindx       #mJy
	    elif t[t['name']==psr[i]]['S1400'][0] != '*':
		s350=float(t[t['name']==psr[i]]['S1400'][0])*(0.35/1.4)**spindx      #mJy
	    else:
		print "cannot calculate s350 for %s" %psr[i]
		continue
	    if t[t['name']==psr[i]]['W50'][0] != '*':
		w50=float(t[t['name']==psr[i]]['W50'][0])     #ms
	    else:
		width_flag=True
	#	print "no catalog width - defaulting to 6\% of period"
		w50=t[t['name']==psr[i]]['P0'][0]*60   #ms
	    dm=float(t[t['name']==psr[i]]['DM'][0])            #pc/cc
	    p0="%.4f" %t[t['name']==psr[i]]['P0'][0]                   #s
	    if t[t['name']==psr[i]]['temp'][0] != '*':
		temp=23+float(t[t['name']==psr[i]]['temp'][0])     #K
	    else:
		print "no temp value for %s" %psr[i]
		continue
	    ts=10.**(-6.46+0.154*math.log10(dm)+1.07*(math.log10(dm))**2-3.86*math.log10(centfreq/1000.))   #ms
	    tdm=8.3*10**6*dm*df/centfreq**3     #ms
	    weff=math.sqrt(w50**2+tdm**2+ts**2) #ms
	    w50b=w50/(20.*float(p0))    #bins
	    weffb=weff/(20.*float(p0))  #bins
    fits=glob('/lustre/cv/projects/GBNCC/amcewen/%s_temp/*%s*fits' %(psr[i],beam[i]))[0]
    poln=float(sproc.Popen(["readfile %s | grep polns | awk '{print $5}'" %fits],stdout=sproc.PIPE,shell=True).communicate()[0][0])
    prep=glob('/lustre/cv/projects/GBNCC/amcewen/%s_temp/prepfold*%s*' %(psr[i],beam[i]))[0] 
    numerator=sproc.Popen(["grep points %s | tail -1" %prep],stdout=sproc.PIPE,shell=True).communicate()[0].split(' ')[6][0:-1]
    try:
	numerator=float(numerator)
    except:
	numerator=float(numerator[1])
    denominator=sproc.Popen(["grep points %s | head -1" %prep],stdout=sproc.PIPE,shell=True).communicate()[0].split(' ')[-1][0:-1]
    try:
	denominator=float(denominator)
    except:
	denominator=float(denominator[1])
    bandwidth=100.*numerator/denominator               #MHz
    sep=float(sproc.Popen(["grep %s /users/amcewen/GBNCC_work/*update*  | grep %s | grep %s | awk '{print $7}'" %(psr[i],p0,beam[i])],stdout=sproc.PIPE,shell=True).communicate()[0][0:-1])   #deg 
    prof_file=glob('/lustre/cv/projects/GBNCC/amcewen/%s_temp/*%s*bestprof' %(psr[i],beam[i]))[0]
    num,val=np.loadtxt(prof_file, dtype='float', unpack=True)
    lim = val < val.mean()
    for j in range(3):
	lim2 = val < val[lim].mean()+2.5*val[lim].std()
	lim = val < val[lim2].mean()+2.5*val[lim2].std()
	if val[lim].std() == val[lim2].std():
	    break
    W = len(val)-len(val[lim])+0.0
    if W == 0.0:    
	W = 1.0    #bins
    elif W > 1:
	high_val = num[val == val.max()][0]
	for i in range(len(lim)):
	    if i != high_val and lim[i] == False and \
	       lim[i-1] == True and lim[(i+1)%len(lim)] == True:
		lim[i] = True
	W = len(val)-len(val[lim])+0.0        #bins
    pfc=2*math.sqrt(poln*120*bandwidth)
    f=math.exp((-(sep/0.3)**2)/1.5)
    snr_data=(s350/temp)*pfc*math.sqrt((P-W)/W)*f
    snr_catW=(s350/temp)*pfc*math.sqrt((P-w50b)/w50b)*f
    try:
	snr_smearW=(s350/temp)*pfc*math.sqrt((P-weffb)/weffb)*f
    except:
	snr_smearW='*'
    string="%s \t %s \t %.4s \t %s \t %.4s \t %s \t %.7s  %.7s \t %s \t %.7s \t %.4s \t %.7s \t %.4s \t %.7s \n" %(psr[i],beam[i],bandwidth,poln,sep,spindx,s350,temp,W,snr_data,w50b,snr_catW,weffb,snr_smearW)
    if width_flag:
	string = string+" !"
    print string
    out.write(string)
print "Note: ! indicates lack of catalog value for W50, so 6% of the period is used. Also, *s in the snr_smearW column indicate cases where the width is smeared out beyond the period width."
