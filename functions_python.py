import numpy as np
import subprocess as sproc
import matplotlib.pyplot as plt
from astropy.table import Table
from glob import glob
import os
import math


def help():				# this help, dummy
    for ln in sproc.Popen("grep def functions_python.py | grep -v communicate",shell=True,stdout=sproc.PIPE).communicate()[0].split('\n'):
	print ln

def tskypy(gb,gl):                      # calculate system temperature at 350MHz for a given sky position; assumes Tsys-Tgal=23K 
    tskypath = "/Users/aemcewen/GBNCC_paper/jail/tsky.ascii"
    tskylist = []
    with open(tskypath) as f:
        for line in f:
            str_idx = 0
            while str_idx < len(line):
                # each temperature occupies space of 5 chars
                temp_string = line[str_idx:str_idx+5]
                try:
                    tskylist.append(float(temp_string))
                except:
                    pass
                str_idx += 5

    # ensure l is in range 0 -> 360
    b = gb
    if gl < 0.:
        l = 360 + gl
    else:
        l = gl

    # convert from l and b to list indices
    j = b + 90.5
    if j > 179:
        j = 179
    nl = l - 0.5
    if l < 0.5:
        nl = 359
    i = float(nl) / 4.
    tsky_haslam = tskylist[180*int(i) + int(j)]
    # scale temperature and add other sources of heat (23K) before returning 
    return tsky_haslam * (350/408.0)**(-2.6)+23

def nor(l1,l2):                         # return boolean array from NOR combination of boolean lists
    if len(l1)==len(l2):
        l=[True for i in range(len(l1))]
        for i in range(len(l1)):
            if l1[i]==False and l2[i]==False:
                l[i]=False
        return np.array(l)
    else:
        print "lists aren't same length"

def snm_err(x,y,p0s,snm,w,dw,lim,verbose=False):      # calculate error on S/N from errors on sigma, width, individual bins, and off-pulse mean
    p0b=len(x)
    dsig=(snm/y[lim].std())*np.sqrt((p0s/120.)*np.sum(y[lim]-y[lim].mean())**2 + (dw/(2*y[lim].std()*(p0b-w)))**2)
    dwer=snm*dw/w
    dpi=(np.sqrt(p0s/120.))*p0b/np.sqrt(w)
    dpo=p0b/(np.sqrt(w*(p0b-w)))
    dsnm=np.sqrt(dpo**2+dsig**2+dwer**2+dpi**2)
    if verbose:
	print "Mean error term = %2.3g" %dsig
	print "Width error term = %2.3g" %dwer
	print "Bin error term = %2.3g" %dpi
	print "Off-pulse error term = %2.3g" %dpo
    return dsnm

def snm(n,v,lim,p,verbose=False):                     # calculate S/N and error on S/N and print DOF correction; period in ms
    power, factor = 1.806, 0.96
    dt_per_bin=p/(float(len(v))*8.192e-05)  # period in ms
    dof_corr=dt_per_bin*factor*(1.0+dt_per_bin**power)**(-1./power)
    #print "DOF correction = %f" %dof_corr
    w=len(lim)-sum(lim)
    snm=dof_corr*sum(v-v[lim].mean())/(v[lim].std()*np.sqrt(w))
    dsnm=snm_err(n,v,p/1000.,snm,w,0.5,lim,verbose=verbose)
    return snm,dsnm,dof_corr

def s350(snm,temp,bw,w,p,theta,dsnm,beta=1.1,dtheta=0.5):   # calculate flux density and error at 350MHz 
    s350v=(temp*beta*snm/(2*np.sqrt(2*120*bw*(p/w-1))))*np.exp((theta/0.3)**2/1.5)
    ds=s350v*np.sqrt((10/temp)**2+(dsnm/snm)**2+(5/bw)**2+(dtheta/(2*w*(1-w/p)))**2) #assumes 10K temperature unc, 5MHz bw unc, 0.5 deg theta unc
    return s350v,ds

def shift(l,n):                         # shift list by n places
    return l[n:]+l[:n]

def plot(n,v,name,surv):   		# plot profiles with name and shifted so that peak is at 1/4 of profile length
    while v.tolist().index(v.max())!=len(n)/4:
	v=v.tolist()
	v=shift(v,1)
	v=np.array(v)
    plt.plot(n,v)
    plt.xticks([])
    plt.yticks([])
    plt.text(int(len(v)/3),int(v.min()+3*(v.max()-v.min())/4),name,fontsize=6)
#    plt.text(int(len(v)/3),int(v.min()+3*(v.max()-v.min())/5),surv,fontsize=6)
#    plt.title(name)


def jnamesort(l):   			# returns list of bestprof files that are sorted by Jname (assuming that files begin with Jname)
    tmp=[]
    for i in l:
	tmp.append([i.split('/')[-1].split('_')[0],i.split('/')[-1].split('_')[1]])
    ls=[]
    tmp.sort()
    for i in range(len(tmp)):
	string=''
	for j in range(len(l[i].split('/'))-1):
	    if j==0:
		string=l[i].split('/')[j]
	    else:
		string=string+'/'+l[i].split('/')[j]
	string=string+'/'
	ls.append(glob(string+"*%s*%s*bestprof" %(tmp[i][0],tmp[i][1]))[0])
    return ls

def catalog(psr,flags):   		# function to call the pulsar catalog from python 
    return sproc.Popen(["/Users/aemcewen/Downloads/psrcat_tar/psrcat %s %s | grep %s" %(flags,psr,psr)],stdout=sproc.PIPE,shell=True).communicate()[0].split()


def wsmr(w,dm,verbose=False):   			# calculate the bin smearing from ts and td using Bhat et. al results; specifically for GBNCC
      cent_freq=0.350 #GHz
      df=0.0244140625 #MHz
      tdm=8.24e-6*float(dm)*df/cent_freq**3  #s 
      ts=10**(-3-6.46+0.154*np.log10(float(dm))+1.07*(np.log10(float(dm)))**2-3.86*np.log10(float(cent_freq))) #s
      if verbose:
           print "DM smearing = %2.3g s" %tdm
      	   print "Scattering time = %2.3g s" %ts
      return np.sqrt(w**2+tdm**2+ts**2)

def freq_adj(s,f,spindx):   		# returns flux at frequency f (in MHz) given flux at 350MHz; spindx<0
      return s*(350./f)**(spindx)

def dot(l1,l2):   			# returns dot product of two boolean lists (boolean AND function)
      l1=np.array(l1).astype(int)
      l2=np.array(l2).astype(int)
      if len(l1)==len(l2):
	  l=l1+l2-1
	  for i in range(len(l)):
	     if l[i]!=1:
		l[i]=0
          return l.astype(bool)
      else:
	print "lists aren't same length, idiot"

def tex_error(val_list,err_list):       # prints out values and errors that have been formatted for latex
    vwerr_list=[]
    maxlen=0
    for v in val_list:
        if len(str(v).split('.')[0])>maxlen:
            maxlen=len(str(v).split('.')[0])
    for v,err in zip(val_list,err_list):
        if v=='*' or err=='*':
            vwerr_list.append('$-$')
            continue
        v=str(v)
        err=float(err)
        errpow=int(("%1.2e" %err).split('e')[-1])
        errstr=''
        for c in str(err):
	    flg=False
            if errpow<0:
                if c=='0' or c=='.':
                    continue
                errstr=errstr+str(c)
                if len(errstr)==2:
                    break
            elif errpow==0:
		if c=='1':
		    if len(errstr)<=1:
		        errstr=errstr+c
		    if len(errstr)==1:
		        continue
		    else:
			break
                if c=='.':
                    continue
                errstr=errstr+c
		break
	    else:
		flg=True
		errstr='6'
		break
        if errstr[0]!='1':
            errstr=errstr[0]
        else:
            errpow=errpow-1
        if flg:
	    vwerr=str(int(np.round(float(v),decimals=-1)))+'('+str(int(np.round(err,decimals=-1)))+')'
	else:
	    vwerr=v.split('.')[0]+'.'+v.split('.')[-1][0:-errpow]+'('+errstr+')'
	if len(vwerr.split('.'))>1:
	    for i in range(maxlen-len(vwerr.split('.')[0])):
	        vwerr='\\phn'+vwerr
	else:
	    for i in range(maxlen-len(vwerr.split('(')[0])):
		vwerr='\\phn'+vwerr
	vwerr_list.append(vwerr)
    return vwerr_list

def inv(l):  				# prints inverse of boolean array (boolean NOT)
    li=np.zeros(len(l)).astype(bool)
    for i in range(len(li)):
	li[i]=not l[i]
    return li

def detectable(p,dm,npol=2,snrmin=5,g=2,tint=120.,dc=0.08,beta=1.1,bw=67):    # roughly determines if period/dm combination is detectable in GBNCC
    i=0
    flg=False
    #tsys=23+34.4*(350./408)**-2.6
    tsys=23+20*(350./408)**-2.75
    pfc=snrmin*beta*tsys/(g*np.sqrt(npol*tint*bw))
    x=np.arange(0.00005,8,0.0001)
    for v in x:
        sensline=pfc*np.sqrt(wsmr(v*dc,dm)/(v-wsmr(v*dc,dm)))
        if str(sensline)!='nan' and not flg:
            print "minimum period and flux = %2.3g s, %2.3g mJy" %(v, sensline)
            flg=True
        if flg: 
            break
    psrsensline=pfc*np.sqrt(wsmr(p*dc,dm)/(p-wsmr(p*dc,dm)))
    if str(psrsensline)!='nan':
        print "for given period and dm (%2.3g s, %2.3g pc/cc), minimum flux is %2.3g mJy" %(p,dm,psrsensline)
        print "smeared pulse width/period (undetectable if >=1) = %2.3g" %(wsmr(0.06*p,dm)/p)
    else:
        print "given pulsar (%2.3g s, %2.3g pc/cc) is not detectable" %(p,dm)
