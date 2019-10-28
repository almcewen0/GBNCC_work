#@ Calculates absolute (extrapolated) pulse widths for a given fraction of the period
#@ Reads in bestprof files (or all bestprofs in current directory)


import numpy as np
import matplotlib.pyplot as plt
from subprocess import call,PIPE,Popen
from glob import glob
from astropy.table import Table
from functions_python import snm as snm_fnc,s350 as s350_fnc

def shift(l,n):
    return l[n:]+l[:n]

def catalog(psr,flags):
    return Popen(["/Users/aemcewen/Downloads/psrcat_tar/psrcat %s %s | grep %s" %(flags,psr,psr)],stdout=PIPE,shell=True).communicate()[0].split()

def remove_bw(bw,x1,x2):
    res1=(100./4096)*float(x1)
    res2=(100./4096)*float(x2)
    return bw-float(('%f' %(res2-res1)))

def width_SAD(x,y,frac,lim):
    y=y-y[lim].mean()
    rts=[]
    ym=y.max()
    ip=[list(y).index(ym)]
    xm=x[ip[0]]
    ip2=[]
    while y[np.array(ip).max()+1]>frac*ym:
        ip.append(np.array(ip).max()+1)
    while y[np.array(ip).min()-1]>frac*ym:
        ip.append(np.array(ip).min()-1)
    if y[list(y).index(ym)+len(y)/2]>frac*ym:
        ip2=[list(y).index(ym)+len(y)/2]
        ym2=y[ip2[0]]
        xm2=x[ip2[0]]
        while y[np.array(ip2).max()+1]>frac*ym:
            ip2.append(np.array(ip2).max()+1)
        while y[np.array(ip2).min()-1]>frac*ym:
            ip2.append(np.array(ip2).min()-1)
        ip2=np.array(ip2)
        m=y[ip2.min()]-y[ip2.min()-1]
        b=y[ip2.min()]-m*x[ip2.min()]-ym*frac
        rts.append(np.roots([m,b])[0])  
        m=y[ip2.max()+1]-y[ip2.max()] 
        b=y[ip2.max()]-m*x[ip2.max()]-ym*frac
        rts.append(np.roots([m,b])[0])
    lim2=(np.zeros(len(y))+1).astype(bool)
    ip=np.array(ip)
    ip.sort()
    for i in range(len(y)):
        if i in ip or i in ip2:
            lim2[i]=False
    m=y[ip.min()]-y[ip.min()-1]
    b=y[ip.min()]-m*x[ip.min()]-ym*frac
    rts.append(np.roots([m,b])[0])
    m=y[ip.max()+1]-y[ip.max()]
    b=y[ip.max()]-m*x[ip.max()]-ym*frac
    rts.append(np.roots([m,b])[0])
    w=0.0
    for i in range(len(rts)):
        if i!=0 and (i+1)%2==0:
            w+=rts[i]-rts[i-1]
    return x,y,lim2,w,rts

def width(x,y,frac,lim,manual_bool):
    y=y-y[lim].mean()
    hy=y.max()*frac
    rts=[]
    for i in range(len(y)-1):
        if i>0 and y[i]>hy:
	    if manual_bool and lim[i]==True:
		continue
            if y[i-1]<hy and y[i+1]<hy:
                m=y[i]-y[i-1]
                b=y[i]-m*x[i]-hy
                rts.append(np.roots([m,b])[0])
                m=y[i+1]-y[i]
                b=y[i]-m*x[i]-hy
                rts.append(np.roots([m,b])[0])
            elif y[i-1]<hy:
                m=y[i]-y[i-1]
                b=y[i]-m*x[i]-hy
                rts.append(np.roots([m,b])[0])
            elif y[i+1]<hy:
                m=y[i+1]-y[i]
                b=y[i]-m*x[i]-hy
                rts.append(np.roots([m,b])[0])
    
    w=0.0
    for i in range(len(rts)):
        if i%2==1:
            w+=rts[i]-rts[i-1]
    return x,y,lim,w,rts


cent_freq=350 #MHz
df=0.0244140625 #MHz

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

writefile='psr_widths_120818.txt'
done=np.loadtxt(writefile,usecols=[0],dtype=str)
rfi=np.loadtxt('/Users/aemcewen/GBNCC_paper/chans_rfi.txt',usecols=[0],dtype=str)
profdir='presto_plots/master_112818/good_profs/'
tmp=open(writefile,'a')


frac=np.array([0.1,0.5])
frac_pct=100*frac
for p in resn['name']:
    if p in done:
        continue
    cond=resn['name']==p
    print p
#    snm=resn['snm'][cond]
    temp=resn['temp'][cond]
    bw=float(resn['bw'][cond])
    if p in rfi:
	for chans in Popen('grep %s /Users/aemcewen/GBNCC_paper/chans_rfi.txt' %p, shell=True,stdout=PIPE).communicate()[0].split()[-1].split(','):
	    bw=remove_bw(bw,chans.split(':')[0],chans.split(':')[1])
    p0=resn['p0'][cond][0]
    off=resn['dist'][cond]
    prof=glob('%s%s*bestprof*' %(profdir,p))
    if len(prof)==0:
	print "Missing profile for %s" %p
	print ""
	continue
    prof=prof[0]
    x,y=np.loadtxt(prof,unpack=True)
    y=y-y.mean()
    y=y.tolist()
    while y.index(np.array(y).max())!=len(y)/4:
        y=shift(y,1)
    y=np.array(y)
    lim1=y<frac[0]*y.max()
    lim2=lim1
    x1,y1,lim1,w1,rts1=width_SAD(x,y,frac[0],lim1)
    x2,y2,lim2,w2,rts2=width_SAD(x,y,frac[1],lim2)
    flg=True
    while flg:
	lim2=lim1
        fig,ax=plt.subplots(1,2,figsize=(12,8))
        ax[0].plot(x1,y1,lw=0.6)
        ax[0].axhline(y=frac[0]*y1.max(),label='W'+str(int(frac_pct[0])))
        ax[0].plot(x1[lim1],y1[lim1],lw=0.6)
	ax[0].legend()
	ax[1].plot(x2,y2,lw=0.6)
	ax[1].axhline(y=frac[1]*y2.max(),label='W'+str(int(frac_pct[1])))
	ax[1].plot(x2[lim2],y2[lim2],lw=0.6)
	ax[1].legend()
        fig.suptitle(p) 
        for r in rts1:
            ax[0].axvline(x=r)
	for r in rts2:
	    ax[1].axvline(x=r)
        plt.show()
        resp=raw_input('good width, switch to cutoff mode, shift, or skip (y/n/s/h/k)? ')
#	resp='y'
        if resp=='y':
            w1ms=(p0/len(y))*1000.*w1
            dc1=w1ms/(1000.*p0)
	    w2ms=(p0/len(y))*1000.*w2
	    dc2=w2ms/(1000.*p0)
	    fig,ax=plt.subplots(1,2,figsize=(12,8))
	    ax[0].plot(x1,y1,lw=0.6)
            ax[0].axhline(y=frac[0]*y1.max(),label='W'+str(int(frac_pct[0])))
            ax[0].plot(x1[lim1],y1[lim1],lw=0.6)
	    ax[0].text(0.5*x1.max(),0.5*(y1.max()-y1.min()),'W%i = %2.3f ms' %(frac_pct[0],w1ms))
            ax[0].legend()
            ax[1].plot(x2,y2,lw=0.6)
            ax[1].axhline(y=frac[1]*y2.max(),label='W'+str(int(frac_pct[1])))
            ax[1].plot(x2[lim2],y2[lim2],lw=0.6)
	    ax[1].text(0.5*x2.max(),0.5*(y2.max()-y2.min()),'W%i = %2.3f ms' %(frac_pct[1],w2ms))
            ax[1].legend()
            fig.suptitle(p)
            for r in rts1:
                ax[0].axvline(x=r)
            for r in rts2:
                ax[1].axvline(x=r)
	    plt.savefig('presto_plots/master_112818/width_profs/%s.png' %p)
	    plt.close()
	    snm,dsnm,dof=snm_fnc(x1,y1,lim1,p0*1000.)
	    if off<=0.5:
	        s350,ds=s350_fnc(snm,temp,bw,len(lim1)-sum(lim1),len(lim1),off,dsnm)
	    else:
		print "using 1 deg uncertainty"
		s350,ds=s350_fnc(snm,temp,bw,len(lim1)-sum(lim1),len(lim1),off,dsnm,dtheta=1.0)
            if p0>1.0:
                print "PSR %s\nP0 = %2.3f s\nW%i = %2.3f ms\nDC%i = %2.3f\nW%i = %2.3f ms\nDC%i = %2.3f\nOn-pulse = %i/%i\nS/N = %2.5f(%2.5f)\nDOF adj. = %2.5f\nS350 = %2.5f(%2.5f)" %(p,p0,frac_pct[0],w1ms,frac_pct[0],dc1,frac_pct[1],w2ms,frac_pct[1],dc2,len(lim2)-sum(lim2),len(lim2),snm,dsnm,dof,s350,ds)
		print ""
            else:
                print "PSR %s\nP0 = %2.3f ms\nW%i = %2.3f ms\nDC%i = %2.3f\nW%i = %2.3f ms\nDC%i = %2.3f\nOn-pulse = %i/%i\nS/N = %2.5f(%2.5f)\nDOF adj. = %2.5f\nS350 = %2.5f(%2.5f)" %(p,p0*1000.,frac_pct[0],w1ms,frac_pct[0],dc1,frac_pct[1],w2ms,frac_pct[1],dc2,len(lim2)-sum(lim2),len(lim2),snm,dsnm,dof,s350,ds)
		print ""
	    tmp.write("%s %2.3f %2.3f %2.3f\n" %(p,p0,w1ms,w2ms)) 
	    profabr=prof.split('/')[-1]
	    call("grep '#' %s > ./presto_plots/master_112818/good_profs_withwidth/%s.width" %(prof,profabr),shell=True)
	    tmp2=open("./presto_plots/master_112818/good_profs_withwidth/%s.width" %profabr,'a')
	    for i in range(len(x)):
		tmp2.write("%i %e %s\n" %(x1[i],y1[i],int(lim1[i])))
	    tmp2.close() 
            flg=False
        elif resp=='n':
	    #resp2=='1'
	    #if resp2=='1':
            for i in range(len(lim1)):
                if lim1[i]==False:
                    print i
            while resp!='g':
                resp=raw_input('input comma separated range to flip or g to finish: ')
                if resp!='g':
                    for i in range(int(resp.split(',')[0]),int(resp.split(',')[1])):
                        lim1[i]=not lim1[i]
            x1,y1,lim1,w1,rts1=width(x1,y1,frac[0],lim1,True)
	    x2,y2,lim2,w2,rts2=width(x2,y2,frac[1],lim2,True)
	    #else:
            #    for i in range(len(lim2)):
            #        if lim2[i]==False:
            #            print i 
            #    while resp!='g':
            #        resp=raw_input('input bins to zap or g to finish: ')
            #        if resp!='g':
            #            resp=int(resp)
            #            lim2[resp]=not lim2[resp]
            #    x2,y2,lim2,w2,rts2=width(x2,y2,frac[1],lim2,True) 
	elif resp=='s':
	    resp2=raw_input('which width (1/2/b)?')
	    if resp2=='1':
		x1,y1,lim1,w1,rts1=width(x1,y1,frac[0],lim1,False)
	    elif resp2=='2':
		x2,y2,lim2,w2,rts2=width(x2,y2,frac[1],lim2,False)
	    else:
		x1,y1,lim1,w1,rts1=width(x1,y1,frac[0],lim1,False)
		x2,y2,lim2,w2,rts2=width(x2,y2,frac[1],lim2,False)
	elif resp=='h':
	    shft=int(raw_input("how many places? "))
	    y1=np.array(shift(list(y1),shft))
	    lim1=np.array(shift(list(lim1),shft))
	    y2=np.array(shift(list(y2),shft))
	    lim2=np.array(shift(list(lim2),shft))
	    x1,y1,lim1,w1,rts1=width_SAD(x1,y1,frac[0],lim1)
	    x2,y2,lim2,w2,rts2=width_SAD(x2,y2,frac[1],lim2)
	elif resp=='k':
	    flg=False
	else:
	    pass
