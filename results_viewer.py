import numpy as np
import subprocess as sproc
from glob import glob
import matplotlib.pyplot as plt
from astropy.table import Table,Column


readfile='psrcat_v158_north_v5_zuul_1720pulsars.txt'
#readfile_dm='DM_search_v4_zuul_1720pulsars.txt'

name=sproc.Popen(["grep -v @ %s | grep -v Not | awk '{print $1}' | uniq" %readfile],stdout=sproc.PIPE,shell=True).communicate()[0].split('\n')[0:-1]
#name_dm=sproc.Popen(["grep -v @ %s | grep -v Not | awk '{print $1}' | uniq" %readfile_dm],stdout=sproc.PIPE,shell=True).communicate()[0].split('\n')[0:-1]

tmp=open('gbncc_duty_cycles','w')
tmp.write("#@ Duty cycles for highest S/N detections in GBNCC Data\n")
tmp.write("#@ Name \t P0 \t W50 \t S/N \t Beam \n")

f=[]
psr_p0=[]
t=Table(names=('name','rad','decd','p0','dm','Dist','MJD','BW', 'Spindx', 'Temp','Wdat','Wcat','Wsmr','S/Ne','S/Nc','S/Ns','S/Nm','Flux', 'dFlux', 'Beam'),dtype=('S15','float','float','float','float','float','int','float','S10','S10','float','S10','S10','S10','S10','S10','S10','S10','S10','int'))
#psr_data_dm=[]
#f_dm=[]
#psr_p0_dm=[]
#t_dm=Table(names=('name','rad','decd','p0','dm','dms','Dist','MJD','BW', 'Spindx', 'Temp','Wdat','Wcat','Wsmr','S/Ne','S/Nc','S/Ns','S/Nm','Flux', 'dFlux', 'Beam'),dtype=('S15','float','float','float','float','str','float','int','float','S10','S10','float','S10','S10','S10','S10','S10','S10','S10','S10','int'))


# Finds highest S/N detection
for i in name:
    snrmax=0.0
    for line in sproc.Popen(["grep %s %s | grep -v Not" %(i,readfile)],stdout=sproc.PIPE,shell=True).communicate()[0].split('\n')[0:-1]:
	if float(line.split('\t')[16])>snrmax:
	    snrmax=float(line.split('\t')[16].split()[0]) 
	    x=line.split()
    for i in range(len(x)-1):
	x[i]=x[i].split()[0]
    for i in range(4):
	x.pop(-2)
    t.add_row(np.array(x))


# Width finder
def shift(l,n):
    return l[n:]+l[:n]
cent_freq=350 #MHz
df=0.0244140625 #MHz
prof=[]
width=[]

frac=0.5

for i in range(len(t)):
    x,y=np.loadtxt(glob("/lustre/cv/projects/GBNCC/amcewen/detection_plots/*%s*%s*bestprof" %(t['Beam'][i],t['name'][i]))[0],unpack=True)
    lim=y<y.mean()
    for i in range(3):
	lim2=y<y[lim].mean()+2.5*y[lim].std()
	lim=y<y[lim2].mean()+2.5*y[lim2].std()
	if y[lim].std() == y[lim2].std():
	    break
    y=y-y[lim].mean()
    hy=y.max()*frac
    rts=[]
    if y[0]>hy or y[-1]>hy:
	y=y.tolist()
	y=shift(y,sum(y>hy)+1)
	y=np.array(y)
    for i in range(len(y)-1):
	if i>0 and y[i]>hy:
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
    if len(rts)<2:
	w='*'
	ws='*'
    else:
	w=0.0
	for i in range(len(rts)):
	    if i%2==1:
		w+=rts[i]-rts[i-1]
    if t['p0'][i]!='*':
	ws=w*t['p0'][i]/len(x)
    else:
	ws='*'
    width.append(ws)

t.add_column(Column(width))	









#for i in name_dm:
#    snrmax=0.0
#    for line in sproc.Popen(["grep %s %s | grep -v Not" %(i,readfile_dm)],stdout=sproc.PIPE,shell=True).communicate()[0].split('\n')[0:-1]:
#        dm=float(line.split('\t')[4].split()[0])
#        p0=float(line.split('\t')[3].split()[0])
#        if float(line.split('\t')[17])>snrmax:
#            snrmax=float(line.split('\t')[17].split()[0]) 
#	    fmax=float(line.split('\t')[18].split()[0])
#	    dms=float(line.split('\t')[5].split('(')[0].split()[0])
#	    dms_err=float(line.split('\t')[5].split('(')[1].split(')')[0].split()[0])
#            beammax=line.split('\t')[-1].split()[-1].split()[0#] 
#	psr_p0_dm.append(p0)
#	f_dm.append(fmax)
#        psr_data_dm.append([i, p0, dm, dms, dms_err, snrmax, fmax, beammax])

#snr_ratio=[]
#name_ratio=[]
#beam_ratio=[]
#dmchange_ratio=[]
#for i in name:
#    if i in name_dm:
#	if psr_data_dm[name_dm.index(i)][-1]==psr_data[name.index(i)][-1] and abs(psr_data_dm[name_dm.index(i)][3]-psr_data_dm[name_dm.index(i)][2])<0.15*psr_data_dm[name_dm.index(i)][2]:
#	    snr_ratio.append(psr_data_dm[name_dm.index(i)][5]/psr_data[name.index(i)][3])
#	    name_ratio.append(i)
#	    beam_ratio.append(psr_data_dm[name_dm.index(i)][-1])
#	    dmchange_ratio.append(psr_data_dm[name_dm.index(i)][3]-psr_data_dm[name_dm.index(i)][2])

#print "%s overlapping detections" %len(snr_ratio)
#print "DMS/N / NoDMS/N average = %s" %np.array(snr_ratio).mean()
#print "DMS/N / NoDMS/N standard deviation = %s" %np.array(snr_ratio).std()
#print "Maximum DMS/N / NoDMS/N for %s in beam %s; DM changed by %s " %(name_ratio[snr_ratio.index(np.array(snr_ratio).max())], beam_ratio[snr_ratio.index(np.array(snr_ratio).max())], dmchange_ratio[snr_ratio.index(np.array(snr_ratio).max())])
#print "Minimum DMS/N / NoDMS/N for %s in beam %s" %(name_ratio[snr_ratio.index(np.array(snr_ratio).min())], beam_ratio[snr_ratio.index(np.array(snr_ratio).min())])


#plt.plot(np.log10(psr_p0),np.log10(f),'+',label='NoDMs')
#plt.plot(np.log10(psr_p0_dm),np.log10(f_dm),'+',label='DMs')
#plt.xlabel('P0')
#plt.ylabel('S350')
#plt.legend()
#plt.show()
#resp=raw_input("plot (y/n)?")
#if resp=='y':
#    import matplotlib.pyplot as plt
#    plt.plot(snr_ratio,dmchange_ratio,'+')
#    plt.xlabel('S/N Ratio')
#    plt.ylabel('DM change')
#    plt.show()

