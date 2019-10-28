# (Not-so-)Simple script to find GBNCC beams on zuul and, optionally, process with PRESTO     ###
# Processing requires environment set-up before script execution
# Can make plot of beams (observed, detected, and unobserved) as beamcheck_PSR.png
# Write out information about pulsars found as GBNCC_ATNF_beamcheck_[comp]_[n].txt (or [outrt]_[comp]_[n].txt)
# Combine results from multiple computer systems with combine_beamcheck.py
# Make individual summaries with GBNCC_individual_summaries.py
# Make nice plots of profiles with GBNCC_profile_plots.py
# Written by Renee Spiewak
# Last edit Jun. 5, 2017

# -h,--help     Print help page with flags listed
# --pos         ra and dec in hh:mm:ss.ss (-)dd:mm:ss.ss form
# -p            period in s
# -pd           period derivative in s/s (default: 0)
# --dm          DM in pc/cc
# -f            file name (assumes file has columns of NUMBER, JNAME, RA (deg), DEC (deg), P0 (s), P1 (s/s), DM (cm^-3); include * in unknown ranges)
# --beam        beam number as integer
# -n            number of profile bins (default: 50)
# --nsub        number of channels (default: 128) 
# --npart       number of subintegrations (default: 64)
# --angle       max angular offset between beam and pulsar in deg (default: 0.5 deg)
# --snr         minimum S/N for detection (default: 6) 
# --proc        process files (default: no processing)
# --center      process central (original) beam
# --comp        select computer system (only 'zuul' or 'GB' allowed)
# --pub         Only process published pulsars
# --rfi         Use Scott's options to remove RFI
# -nonew        Do not copy new fits files into directories 
# -o	        Choose output location (default is /lustre/cv/projects/GBNCC/amcewen)
# --par         Provide par file or directory in which par file is stored (if the latter, par file is assumed to be [pulsarJname].par 
# --rew	        Reprocess, i.e. delete previous output
# --outrt       Set root for output file (default is GBNCC_ATNF_beamcheck)
# --psearch     Search period with prepfold
# --dmsearch    Search DMs with prepfold
# --ignorechans Provide file of channels to block for specific pulsars (assumes file has columns of JNAME, CHANS)

#########################################################

import numpy as np
import sys, os
from glob import glob
from subprocess import call,Popen,PIPE
from time import strftime
import argparse as ap
from astropy.coordinates import SkyCoord,Latitude,Longitude 
from astropy.table import Table,Column
from astropy import units as u, constants as c
#import matplotlib.pyplot as plt

# Function to find distance 
def ang_offset(lon1,lat1,lon2,lat2):
	x="%f %f" %(lon1,lat1)
	y="%f %f" %(lon2,lat2)
	return float(SkyCoord(x,unit=("deg","deg")).separation(SkyCoord(y,unit=("deg","deg")))/u.deg)

# Function to cycle items in list (for calculating W50)
def shift(l,n):
    return l[n:]+l[:n]

def proc_args():
    """
    Intelligently read in arguments
    Check for errors

    Returns
    ---------
    dictionary containing relevant values

    """
    pars = ap.ArgumentParser(description="Locate and process data \
on known pulsars")
    pars.add_argument("--comp",action="store",choices=["zuul","GB"],
                      required=True,help="select computer system")
    pars.add_argument("--pos",action="store",
                      help="ra and dec ('hms dms' or degrees)")
    pars.add_argument("-p","--period",action="store",default=1.0,type=float,
                      help="period in s")
    pars.add_argument("-pd","--pdot",action="store",default=0,type=float,
                      help="period derivative in s/s")
    pars.add_argument("--dm",action="store",default=50,type=float,
                      help="DM in pc/cc")
    pars.add_argument("-f","--file",action="store",
                      help="file name (assumes file has columns of \
NUMBER, JNAME, RA, DEC, P0, P1, DM)")
    pars.add_argument("--beam",action="store",
                      help="beam number as integer")
    pars.add_argument("-n","--nbin",action="store",default=200,type=int,
                      help="number of profile bins")
    pars.add_argument("--nsub",action="store",default=128,type=int,
                      help="number of channels")
    pars.add_argument("--npart",action="store",default=64,type=int,
                      help="number of subintegrations")
    pars.add_argument("--angle",action="store",default=0.5,type=float,
                      help="max offset between beam and pulsar in deg")
    pars.add_argument("--snr",action="store",default=6,type=float,
                      help="minimum S/N for detection")
    pars.add_argument("--proc",action="store_true",
                      help="process files")
    pars.add_argument("--center",action="store_true",dest="find_psr",
                      help="process only the closest beam")
    pars.add_argument("--pub",action="store_true",
                      help="only process published pulsars")
    pars.add_argument("--rfi",action="store_true",
                      help="use Scott's options to remove RFI")
    pars.add_argument("--nonew",action="store_true",dest="nonew",help="do not search for .fits files")
    pars.add_argument("-o","--out",action="store",help="output location (default is /lustre/cv/projects/GBNCC/amcewen/)")
    pars.add_argument("--par",action="store",help="parameter file or directory in which parameter file is stored")
    pars.add_argument("--rew",action="store_true",dest="rewrite",help="delete old output and rewrite")
    pars.add_argument("--outrt",action="store",dest="outputroot",help="root for name of summary file")
    pars.add_argument("--dmsearch",action="store_true",dest="dmsearch",help="search DMs with prepfold")
    pars.add_argument("--psearch",action="store_true",dest="psearch",help="search period with prepfold")
    pars.add_argument("--ignorechans",action="store",dest="igchn",help="file with channels to ignore")
    args = vars(pars.parse_args())

    if args['pos'] != None and len(args["pos"].split()) > 1:
        ra_astro = Longitude(args["pos"].split()[0],unit="hourangle")
        dec_astro = Latitude(args["pos"].split()[1],unit="deg")
        args["ra_psr"] = ra_astro.value
        args["dec_psr"] = dec_astro.value

    return args

# Pulsar class to store input information and find other relevant catalog data
class Pulsar:
	def __init__(self,(name,ra,dec,p0,p1,dm)):
		self.name=name
		try:
			ra.split(':')[1]
			self.ra=SkyCoord("%s %s" %(ra,dec),unit=('hourangle','deg')).ra.value
			self.dec=SkyCoord("%s %s" %(ra,dec),unit=('hourangle','deg')).dec.value
		except:
			self.ra=SkyCoord("%s %s" %(ra,dec),unit=('deg','deg')).ra.value
			self.dec=SkyCoord("%s %s" %(ra,dec),unit=('deg','deg')).dec.value
		self.p0=float(p0)
		self.p1=p1
		self.dm=dm
		self.newdm=''
		self.beams=[]
		self.par=''
		self.w50flag=False
		self.params=name,ra,dec,p0,p1,dm	
		self.pos="%s %s" %(ra,dec)
		self.specflgs=''
		if name in specflgpsrs:
		    flg=''
		    for f in Popen("grep %s /users/amcewen/special_flags.txt" %name,shell=True,stdout=PIPE).communicate()[0].split()[1:]:
			flg=flg+' '+f
		    self.specflgs=flg
		if len(cat_dat[cat_dat['name']==name]) != 0:
		    if cat_dat[cat_dat['name']==name]['temp'] != '*':
			self.temp=23+float(cat_dat[cat_dat['name']==name]['temp'])
		    else:
			self.temp=cat_dat[cat_dat['name']==name]['temp']
		    self.w50=cat_dat[cat_dat['name']==name]['W50']
		    if self.w50 == '*':
			self.w50=self.p0*47.8534 
                        self.w50flag=True
		    self.w50b=float(self.w50)/(20.*self.p0)
		    if dm != '*':
			tdm=8.3*10**6*float(dm)*df/cent_freq**3
			ts=10**(-6.46+0.154*np.log10(float(dm))+1.07*(np.log10(float(dm)))**2-3.86*np.log10(float(cent_freq)/1000.))
			self.wsmear=np.sqrt(float(self.w50)**2+ts**2+tdm**2)
			self.wsmearb=float(self.wsmear)/(5.*self.p0)
		    else:
			self.wsmear='*'
			self.wsmearb='*'
		    self.s400=cat_dat[cat_dat['name']==name]['S400'][0]
		    self.s1400=cat_dat[cat_dat['name']==name]['S1400'][0]
		    self.spindx=cat_dat[cat_dat['name']==name]['SPINDX'][0]
                    if self.spindx != '*':
                        if self.s400 != '*':
                            self.s350 = float(self.s400)*(0.35/0.4)**float(self.spindx)
                        elif self.s1400 != '*':
                            self.s350 = float(self.s1400)*(0.35/1.4)**float(self.spindx)
		        else:
			    self.s350 = '*'
                    elif self.s400 != '*' and self.s1400 != '*':
                        self.s350=float(self.s400)**(0.35/0.4)**(np.log(float(self.s400)/float(self.s1400))/-1.253)
			self.spindx=np.log(float(self.s400)/float(self.s1400))/-1.253
                    elif self.s400 != '*':
                        self.s350=float(self.s400)*((0.35/0.4)**(-1.7))
                    elif self.s1400 != '*':
                        self.s350=float(self.s1400)*((0.35/1.4)**(-1.7))
                    else:
                        self.s350='*'
		else:
		    self.temp='*'
		    self.w50=self.p0*47.8534       # calculating width based on average value of W50/p0 from catalog
		    self.w50flag=True
		    self.w50b=float(self.w50)/(20.*self.p0)
                    tdm=8.3*10**6*float(dm)*df/cent_freq**3
                    ts=10**(-6.46+0.154*np.log10(float(dm))+1.07*(np.log10(float(dm)))**2-3.86*np.log10(float(cent_freq)/1000.))
                    self.wsmear=np.sqrt(float(self.w50)**2+ts**2+tdm**2)
                    self.wsmearb=float(self.wsmear)/(5.*self.p0)
		    self.s350='*'
		    self.s400='*'
		    self.s1400='*'
		    self.spindx='*'
		if self.dec > -40 - ang_max:
			self.north=True
		else:
			self.north=False 

	def change_bins(self,nbins):
		if self.w50 != '*':
		    self.w50b=float(self.w50b)*nbins/200.
		    self.wsmearb=float(self.wsmearb)*nbins/200.		
	def beam_numbers(self):
		Str=''
		for i in self.beams:
			Str=Str+"%s " %i.num
		print Str
	def set_pos(self,(ra,dec)):
		self.ra=ra
		self.dec=dec

# Beam class to store GBNCC pointing data
class Beam:
	def __init__(self,beam):
		self.name=beam
		self.num=beam.strip('GBNCC')
		self.mask=[]
                self.maskfrac=[]
		self.fits=[]
		self.mjd=[]
		self.ra=pointings[pointings['pointing']==beam]['RAdeg'][0]
		self.dec=pointings[pointings['pointing']==beam]['Decdeg'][0]
		self.pos="%s %s" %(self.ra,self.dec)
	def add_mask(self,mask):
		self.mask.append(mask)	
                if self.name in mask_fractions['pointing']:
                    self.maskfrac.append(mask_fractions[mask_fractions['pointing']==self.name]['mask_frac'][0])
	def add_fits(self,fits):
		self.fits.append(fits)
		self.mjd.append(int(fits.split('_')[1]))
	def ang_off(self,ra,dec):
		self.off=ang_offset(self.ra,self.dec,ra,dec)

# Begin full program 
args = proc_args()

p0_psr = args["period"]
pd_psr = args["pdot"]
dm_psr = args["dm"]
file_psr = args["file"]
beam_psr = args["beam"]
nbin = args["nbin"]
nsub = args["nsub"]
npart = args["npart"]
ang_max = args["angle"]
snr_min = args["snr"]
proc_psr = args["proc"]
find_psr = args["find_psr"]
use_comp = args["comp"]
f_pub = args["pub"]
rfi_fil = args["rfi"]
par_file = args["par"]
rewrite = args["rewrite"]
outrt = args["outputroot"]
dmsearch = args["dmsearch"]
psearch = args["psearch"]
nonew = args["nonew"]
igchn = args["igchn"]

# Create name of pulsar from command line input
if "ra_psr" in args:
    ra_psr = args["ra_psr"]
    dec_psr = args["dec_psr"]
    psr=str(SkyCoord(ra_psr,dec_psr,unit=('hourangle','deg')).to_string('hmsdms').strip('s'))
    psr_name = "J%s%s%s%s" %(psr.split()[0].split('h')[0], psr.split()[0].split('h')[1].split('m')[0],\
	 psr.split()[1].split('d')[0], psr.split()[1].split('d')[1].split('m')[0])

# Running code with --dmsearch and not --rew will pull "new DM" from old data that did not neccesarily search DMs
if (dmsearch or psearch) and not rewrite:
    print ""
    print "WARNING: Searching is on, but rewrite is not. Previously calculated data will not be recreated with search flag(s)."
    resp = raw_input("Change rewrite flag (--rew) to True (y/n)?\n")
    if resp == 'y':
	args["rewrite"] = True
	rewrite = args["rewrite"]

# Survey specs and directory structure
renee_dir = "/users/rspiewak/pulsars/"
cent_freq=350. #MHz
df=0.0244140625  #MHz
yflg=False

tskypath = glob('/users/amcewen/GBNCC_work/tsky.ascii')[0]	# method for calculating Haslam sky temperature
tskylist=[]
with open(tskypath) as f:
    for line in f:
	str_idx=0
	while str_idx<len(line):
	    temp_string=line[str_idx:str_idx+5]
	    try:
		tskylist.append(float(temp_string))
	    except:
		pass
	    str_idx+=5

def tskypy(ra,dec):
    b=SkyCoord(ra,dec,unit=('deg','deg')).transform_to('galactic').b.value
    l=SkyCoord(ra,dec,unit=('deg','deg')).transform_to('galactic').l.value
    if l<0:
	l=360.+l
    j=b+90.5
    if j>179:
	j=179
    nl=l-0.5
    if l<0.5:
	nl=359
    i=float(nl)/4.
    tsky_haslam=tskylist[180*int(i)+int(j)]
    return tsky_haslam*(350/408.0)**(-2.6)

if "out" in args:
    work_dir = args["out"]
    if work_dir[-1]!='/':
	work_dir=work_dir+'/'
else:
    work_dir = '/lustre/cv/projects/GBNCC/amcewen/'
if use_comp == "zuul":
    data_dir = "/lustre/cv/projects/GBNCC/"
elif use_comp == "GB":
    data_dir = "/lustre/pulsar/survey/AGBT09C_057/"

# Get info on all GBNCC past and future beams from Scott Ransom's file
if not os.path.isfile('./GBNCC_pointings.fits'):
	filename='/lustre/cv/projects/GBNCC/amcewen/GBNCC_posns_by_dec_ALLGBTSKY.txt'
	np.fromfile(filename)
	t=Table.read(filename,format='ascii',names=('pointing','RA','Dec'))
	s=SkyCoord(pointings['RA'],pointings['Dec'],unit=('deg','deg'))
	t.add_column(Column(s.ra,name='RAdeg'))
	t.add_column(Column(s.dec,name='Decdeg'))
	t.write('GBNCC_pointings.fits')
pointings=Table.read('GBNCC_pointings.fits')
mask_fractions=Table.read('/users/amcewen/GBNCC_work/GBNCC_mask_fraction.fits')
cat_dat=Table.read('/users/amcewen/GBNCC_work/psrcat_v159_full.fits')

beam_rejects = np.array([20051,19907,10451,19906,20312,13691,17237,19064,72172,80142,83626,114632,115242,120582])
specflgpsrs=np.loadtxt('/users/amcewen/special_flags.txt',usecols=[0],dtype=str)

print ""
# Check file information (assumes order of data in file)
if not file_psr is None:
    p_coord=[]
    if file_psr.split('.')[-1]!='fits':
    	psr_tbl=Table.read(file_psr,format='ascii')
    else:
	psr_tbl=Table.read(file_psr)
    y=[]
    for dec,p0,dm in zip(psr_tbl[psr_tbl.colnames[3]],psr_tbl[psr_tbl.colnames[4]],psr_tbl[psr_tbl.colnames[6]]):
	if dec>-40 and p0!='*' and dm!='*':
	    y.append(True)
	else:
	    y.append(False)
    print "%d/%d valid pulsars in survey area. Finding closest beams..." %(sum(y),len(psr_tbl))
    if sum(y)==0:
	exit()
# If no file is given, read pulsar from input parameters
else:
    if p0_psr!='*' and dm_psr!='*':
        pulsar=Pulsar((psr_name,ra_psr,dec_psr,p0_psr,pd_psr,dm_psr))
    if pulsar.north:
	print "%s in survey area. Finding closest beams..." %psr_name
	psr_tbl=Table(rows=[[psr_name,ra_psr,dec_psr,p0_psr,pd_psr,dm_psr]])
	yflg=True
	y=np.array([0,1])
    else:
	print "%s not valid and/or in survey area." %psr_name
	exit()

# Create cumulative output file
if "outputroot" in args:
    summary_file = open("%s%s_%s_%dpulsars.txt" %(work_dir,outrt,use_comp,sum(y)), 'w')
else:
    summary_file = open("%sGBNCC_ATNF_beamcheck_%s_%dpulsars.txt" %(work_dir,use_comp,sum(y)), 'w')
summary_file.write("#@ Using S/N > %s for detection threshold on %s\n" %(snr_min,use_comp))
if str(beam_psr) == "None" and not file_psr is None:
    summary_file.write("#@ Using beams within %.1f deg of pulsars in %s \n" %(ang_max,file_psr))
elif str(beam_psr) != "None" and not file_psr is None:
    summary_file.write("#@ Using beam #%s to search for pulsars from %s \n" %(beam_psr,file_psr))
elif str(beam_psr) != "None" and file_psr is None:
    summary_file.write("#@ Using beam #%s to search for pulsar %s \n" %(beam_psr,pulsar.name))
elif str(beam_psr) == "None" and file_psr is None:
    summary_file.write("#@ Using beams within %.1f deg of pulsar %s \n" %(ang_max,pulsar.name))
if 'pub' in args:
    summary_file.write('#@ Only searched for published pulsars\n')
summary_file.write("#@ Flux adjusted for offset from beam center\n")
summary_file.write("#@ Exclamation points next to Wcat and Wsmr indicate cases where W50 is not published, so the average ratio of W50/P0 is used\n")
summary_file.write("#@ S/Ne = expected S/N estimated from published flux values using width calculated from data\n")
summary_file.write("#@ S/Nc = expected S/N calculated using catalog W50\n") 
summary_file.write("#@ S/Ns = S/N calulated using Weff (smeared out)\n")
summary_file.write("#@ S/Nm = measured S/N \n")
summary_file.write("#@ Flux estimated from measured S/N, uncertainty estimated from dependencies\n")
summary_file.write("#@ File generated %s (ET)\n" %strftime("%Y-%m-%d %H:%M:%S"))
if dmsearch:
    summary_file.write('#@ DM Searching - DMs is resulting DM from search\n')
    summary_file.write("#@  PSRJ  \t RA \t Dec \t P0 \t\t DM \t DMs \t\t Dist  \t MJD     BW\t Spindx \
Temp \t Wdat  \t\t Wcat \t Wsmr \t S/Ne \t S/Nc \tS/Ns\tS/Nm \t   Flux\t\tdFlux\t\tStatus \n")
    summary_file.write("#@\t\t(deg)\t (deg) \t (s) \t\t (p/cc)\t (p/cc)\t\t (deg)\t (day)\t(MHz) \t\t (K) \
\t (ms)  \t\t (ms)  \t  (ms)   \t\t\t\t   (mJy)\t(mJy)  \n")
else:
    summary_file.write("#@  PSRJ  \t RA \t Dec \t P0 \t\t DM \t Dist  \t MJD     BW\t Spindx  \
Temp \tWdat  \tWcat \tWsmr \t S/Ne \t S/Nc \t S/Ns\tS/Nm \t   Flux\t\t dFlux\t\tStatus \n")
    summary_file.write("#@\t\t(deg)\t (deg) \t (s) \t\t (p/cc)\t (deg)\t (day)\t(MHz)\t\t (K) \
\t(s)    \t(s)    \t(s)    \t\t\t\t\t   (mJy)\t (mJy)  \n")

# Tracking number for final output
name_detect = []
num_detect = 0
num_unobs = 0

print "Checking S/N for profiles -- threshold = %.1f" %snr_min

if yflg:
    y=[0]
# Find observations and, optionally, process with PRESTO
for i in range(len(psr_tbl[y])):
    if yflg:
	psr=pulsar
    else:
        psr=Pulsar((psr_tbl[y][i][1],psr_tbl[y][i][2],psr_tbl[y][i][3],psr_tbl[y][i][4],psr_tbl[y][i][5],psr_tbl[y][i][6]))
    p_coord=SkyCoord(psr.ra,psr.dec,unit=('deg','deg'))
    if str(beam_psr) == "None":
	if not find_psr:
            for j in range(len(pointings[p_coord.separation(SkyCoord(pointings['RAdeg'],pointings['Decdeg']))<ang_max*u.deg]['pointing'])):
                if not int(pointings[p_coord.separation(SkyCoord(pointings['RAdeg'],pointings['Decdeg']))<ang_max*u.deg]['pointing'][j].strip('GBNCC')) in beam_rejects:
                    psr.beams.append(Beam(pointings[p_coord.separation(SkyCoord(pointings['RAdeg'],pointings['Decdeg']))<ang_max*u.deg]['pointing'][j]))
	else: 
	    psr.beams.append(Beam(pointings[p_coord.separation(SkyCoord(pointings['RAdeg'],pointings['Decdeg']))==p_coord.separation(SkyCoord(pointings['RAdeg'],pointings['Decdeg'])).min()][0][0]))
    elif os.path.isfile(beam_psr):
	psr.beams.append(Beam('GBNCC'+Popen("grep %s %s | awk '{print $2}'" %(psr.name,beam_psr),shell=True,stdout=PIPE).communicate()[0].split('\n')[0]))
    elif int(beam_psr) not in beam_rejects:
	beam_psr='GBNCC'+beam_psr
	psr.beams.append(Beam(beam_psr))
    print ""
    print ""
    print psr.name
    psr.beam_numbers()
    if psr.name[0] == "C": 
	unpub_psr = True
    else:
	unpub_psr = False
    if f_pub and unpub_psr:
	continue     # Skip unpublished pulsars when flag used
    if psr.temp == '*':
	psr.temp=tskypy(psr.ra,psr.dec)     #Calculate otherwise unknown sky temperature from Haslam paper
    if str(par_file) == "None":   # Unless a .par file/directory is given, search through repositories for match
	if len(glob('%sparfiles/%s_short.par' %(renee_dir,psr.name)))>0:
	    psr.par=glob('%sparfiles/%s_short.par' %(renee_dir,psr.name))[0]
	if len(glob(data_dir+'amcewen/pars/%s*' %psr.name))>0:
	    psr.par=glob('/lustre/cv/projects/GBNCC/amcewen/pars/%s*' %psr.name)[0]
	if psr.par == '' and len(psr.name)>=10:
	    print "Cannot find par file in %sparfiles/ or in %samcewen/pars/" %(renee_dir,data_dir)
    else:   # Checking given .par file/directory
	if os.path.isdir(par_file):
    	    try:
		psr.par=glob(par_file+"/*%s*par" %psr.name)[0]
	    except:
		print "Can't find par file for %s in %s" %(psr.name,par_file)
	else:
	    if os.path.isfile(par_file):
		psr.par=par_file
	    else:
		print "Invalid par file"
    os.chdir(work_dir)
    if not os.path.isdir("%s_temp" %psr.name):      # Creating directories for pulsars if they don't exist (and replacing them with empty directories if --rew is set)
	os.mkdir("%s_temp" %psr.name)
    elif not os.path.isdir("%s_temp" %psr.name):
	print "Skipping new pulsar: "+psr.name
	continue
    elif rewrite:
	call("rm -r %s%s_temp" %(work_dir,psr.name),shell=True)
	call("rm %sdetection_plots/*" %work_dir,shell=True)
	os.mkdir("%s%s_temp" %(work_dir,psr.name))
    os.chdir("%s_temp/" %psr.name)
    if not os.path.isfile('%s_short.par' %psr.name) and not psr.par == '':      # Linking .par file to directory
	call('cp -sf %s .' %psr.par,shell=True)
    for beam_cand in psr.beams:
	if os.getcwd().split('/')[-1] != "%s_temp" %psr.name:
	    os.chdir(work_dir+"%s_temp/" %psr.name)
	beam_cand.ang_off(psr.ra,psr.dec)              # Calculate the angular offset between pulsar and beam and store this value in the beam metadata
	if len(glob('*%s*fits' %beam_cand.name))==0:
	    if not nonew:
		fitsfound=False
		dirs=[data_dir,'/lustre/lard/projects/GBNCC/']
		j=0
		for fits_list in [glob('%s2*/gu*G????%s_*fits' %(dirs[0],beam_cand.num)),glob('%s2*/gu*G????%s_*fits' %(dirs[1],beam_cand.num))]:
	        #fits_list=glob('%s2*/gu*G????%s_*fits' %(data_dir,beam_cand.num))  # Search for GBNCC .fits files 
	            if len(fits_list)==1:
	                call("ln -sf %s ." %fits_list[0],shell=True)
		        beam_cand.add_fits(fits_list[0].split('/')[-1])
			fitsfound=True       
	            elif len(fits_list)>1:
		        print "Multiple files for pointing #%s" %beam_cand.num
		        for i in fits_list:
	                    call("ln -sf %s ." %i,shell=True)
		            beam_cand.add_fits(i.split('/')[-1])
			fitsfound=True
                    else:
			fitsdir=''
                        print "Cannot find beam #%s in %s for PSR %s" %(beam_cand.num, dirs[j], psr.name)
		    j+=1
		if not fitsfound:
		    if dmsearch:
                    	summary_file.write("%-12s\t %-5.5s \t %-5.5s \t %-8.8s \t %-.5s \t * \t\t %-.5s \t * \t *  \t %-.4s \t %-.5s \t * \t\t * \%-.3s \t %-.3s \t * \t * \t  * \t   * \t\t *  \t\tBeam %6s Not Found\n" \
                        %(psr.name,psr.ra,psr.dec,psr.p0,psr.dm,beam_cand.off,psr.spindx,psr.temp,psr.w50,psr.wsmear,beam_cand.num))
		    else:
		        summary_file.write("%-12s\t %-5.5s \t %-5.5s \t %-8.8s \t %-.5s \t %-.5s \t * \t * \t %-.4s \t %-.5s \t * \t\t * \%-.3s \t %-.3s \t * \t * \t  * \t   * \t\t *  \t\tBeam %6s Not Found\n" \
                        %(psr.name,psr.ra,psr.dec,psr.p0,psr.dm,beam_cand.off,psr.spindx,psr.temp,psr.w50,psr.wsmear,beam_cand.num))
                    num_unobs += 1
                    os.chdir(work_dir)
                    continue
	    else:
		print "Cannot find beam #%s on %s for PSR %s" %(beam_cand.num, use_comp, psr.name)
		if dmsearch:
		    summary_file.write("%-12s\t %-5.5s \t %-5.5s \t %-8.8s \t %-.5s \t * \t\t %-.5s \t * \t * \t %-.4s \t %-.5s \t * \t\t * \t \
%-.3s \t %-.3s \t * \t * \t  * \t   * \t\t *  \t\tBeam %6s Not Found\n" \
		    %(psr.name,psr.ra,psr.dec,psr.p0,psr.dm,beam_cand.off,psr.spindx,psr.temp,psr.w50,psr.wsmear,beam_cand.num))
		else:
		    summary_file.write("%-12s\t %-5.5s \t %-5.5s \t %-8.8s \t %-.5s \t %-.5s \t * \t * \t %-.4s \t %-.5s \t * \t\t * \t \
%-.3s \t %-.3s \t * \t * \t  * \t   * \t\t * \t\tBeam %6s Not Found\n" \
		    %(psr.name,psr.ra,psr.dec,psr.p0,psr.dm,beam_cand.off,psr.spindx,psr.temp,psr.w50,psr.wsmear,beam_cand.num))
		num_unobs += 1
		os.chdir(work_dir)
		continue
	else:
	    for i in glob('guppi*%s*fits' %beam_cand.name):
	        beam_cand.add_fits(i)
	# Process files to find detections
	rfi_std = "-time 1 -timesig 2 -freqsig 3"
	for fits,mjd in zip(beam_cand.fits,beam_cand.mjd):
	    if os.getcwd().split('/')[-1] != "%s_temp" %psr.name:
		os.chdir("%s_temp/" %psr.name)
	    if "2bit" in fits and rfi_fil == True:		# Setting various flags for rfifind and prepfold
	        raw_opt = " -noscales -noweights"
	        rfi_std = "-time 1 -freqsig 3"
	    else:
	        raw_opt = " "
	    if mjd > 56710 and rfi_fil == True:
	        rfi_opt = " -zapchan 2470:3270"
	    else:
	        rfi_opt = " "
	    while len(glob('%s%s_temp/*%s*%s*rfifind.mask' %(work_dir,psr.name,mjd,beam_cand.num)))==0 \
	        or len(glob('%s%s_temp/*%s*%s*rfifind.stats' %(work_dir,psr.name,mjd,beam_cand.num)))==0:               # If current directory doesn't have mask/stats files, check mask repository
	        if len(glob(data_dir+'amcewen/mask_files/*/*%s*%s*rfifind.mask' %(mjd,beam_cand.num)))!=0 \
		    and len(glob(data_dir+'amcewen/mask_files/*/*%s*%s*rfifind.stats' %(mjd,beam_cand.num)))!=0:
		    mask_dir=glob(data_dir+'amcewen/mask_files/*/*%s*%s_*rfifind.mask' %(mjd,beam_cand.num))[0]
		    stats_dir=glob(data_dir+'amcewen/mask_files/*/*%s*%s_*rfifind.stats' %(mjd,beam_cand.num))[0]
		    call('ln -sf %s .' %mask_dir,shell=True)
		    call('ln -sf %s .' %stats_dir,shell=True)
		if len(glob('%s%s_temp/*%s*%s*rfifind.mask' %(work_dir,psr.name,mjd,beam_cand.num)))!=0:
		    break
		lowlim=int(np.floor(float(beam_cand.num)/5000)*5000)
		highlim=int(np.ceil(float(beam_cand.num)/5000)*5000)
		msk_tar=glob(data_dir+'amcewen/mask_files/GBNCC_%s-%s*tar.gz' %(lowlim,highlim))[0]		    # Next, check my archives
		os.chdir(data_dir+'amcewen/mask_files/')
		print "attempting to untar %s mask/stats files from tarball in %samcewen/mask_files/" %(beam_cand.name,data_dir)
		try:
		    call('tar -xvf %s \*GBNCC%s_\*' %(msk_tar,beam_cand.num),shell=True)
		except:
		    print "failed to untar mask file"
		if len(glob(data_dir+'amcewen/mask_files/*/*%s*%s_*.mask' %(mjd,beam_cand.num)))>0 and len(glob(data_dir+'amcewen/mask_files/*/*%s*%s_*.stats' %(mjd,beam_cand.num)))>0:
		    call('ln -sf %samcewen/mask_files/*/*%s*%s_*.mask %s%s_temp/' %(data_dir,mjd,beam_cand.num,work_dir,psr.name),shell=True)
		    call('ln -sf %samcewen/mask_files/*/*%s*%s_*.stats %s%s_temp/' %(data_dir,mjd,beam_cand.num,work_dir,psr.name),shell=True)
		if len(glob('%s%s_temp/*%s*%s*rfifind.mask' %(work_dir,psr.name,mjd,beam_cand.num)))!=0:
		    break
		try:
		    mask_dir=glob(data_dir+'20*/*%s*GBNCC%s*' %(mjd,beam_cand.num))[0].split('/')[5]        # If they aren't found there, try to get them from Scott's archives   
		except:
		    mask_dir=glob('/lustre/lard/projects/GBNCC/20*/*%s*GBNCC%s*' %(mjd,beam_cand.num))[0].split('/')[5]
		if len(glob(data_dir+mask_dir+'/*rfi*tar*'))>0:
		    print "untarring mask file"
		    os.chdir(data_dir+'amcewen/mask_files/')
		    for tar_file in glob(data_dir+mask_dir+'/*rfi*tar*'):
			if len(glob('*%s*%s*rfifind.mask' %(mjd,beam_cand.num)))==0 or len(glob('*%s*%s*rfifind.mask' %(mjd,beam_cand.num)))==0:
			    call('tar -xf %s \*GBNCC%s_\*find.*' %(tar_file,beam_cand.num), shell=True)
			    call('ln -sf *%s*GBNCC%s*rfifind.mask %s*%s*/' %(mjd,beam_cand.num,work_dir,psr.name),shell=True)
			    call('ln -sf *%s*GBNCC%s*rfifind.stats %s*%s*/' %(mjd,beam_cand.num,work_dir,psr.name),shell=True)
		os.chdir('%s%s_temp/' %(work_dir,psr.name))
		print "Running rfifind for %s beam candidate %s..." %(psr.name,beam_cand.num)
		rfi_out = open('rfifind_%s_%s_output.txt' %(mjd,beam_cand.num),'w')
		rfi_out.write("nice rfifind %s%s%s -o %s_%s_%s %s"
			%(rfi_std,raw_opt,rfi_opt,psr.name,mjd,beam_cand.num,fits))
		p = Popen("nice rfifind %s%s%s -o %s_%s_%s %s"
			%(rfi_std,raw_opt,rfi_opt,psr.name,mjd,beam_cand.num,fits),
		        stdout=rfi_out, shell=True)
		p.wait()
		rfi_out.flush()
		rfi_out.close()
		del p
		if not os.path.isdir('%samcewen/mask_files/GBNCC_%s-%s_mask' %(data_dir,lowlim,highlim)):
		    call('mkdir %samcewen/mask_files/GBNCC_%s-%s_mask' %(data_dir,lowlim,highlim),shell=True)
		call("mv *%s*%s_*d.mask %samcewen/mask_files/GBNCC_%s-%s_mask" %(mjd,beam_cand.num,data_dir,lowlim,highlim),shell=True)
                call("mv *%s*%s_*d.stats %samcewen/mask_files/GBNCC_%s-%s_mask" %(mjd,beam_cand.num,data_dir,lowlim,highlim),shell=True)
		call("mv *%s*%s_*rfifind.* %samcewen/mask_files/not_mask_stats"  %(mjd,beam_cand.num,data_dir),shell=True)
		call('ln -sf %samcewen/mask_files/*/*%s*%s_*d.mask .' %(data_dir,mjd,beam_cand.num), shell=True)        # Move mask/stats files to the mask repository and link to them
		call('ln -sf %samcewen/mask_files/*/*%s*%s_*d.stats .' %(data_dir,mjd,beam_cand.num), shell=True) 

	    os.chdir(work_dir+'%s_temp' %psr.name)
	    beam_cand.add_mask(glob('*%s*%s*rfifind.mask' %(mjd,beam_cand.num))[0])
	    if len(beam_cand.fits) == len(beam_cand.maskfrac) and beam_cand.maskfrac[-1] != '':                 # Calculate mask fraction (if it isn't already reported in Will's file)
	        bandwidth=100.*(1-beam_cand.maskfrac[-1])
	    if not dmsearch:				 # Setting various prepfold flags 
	        flag_search = '-nopsearch -nodmsearch'
	        if (len(psr.name) < 10 and psr.p0 < 0.1) or psearch:
		    flag_search = '-nodmsearch'
	    elif (len(psr.name) < 10 and psr.p0 < 0.1) or psearch:
	        flag_search=""
	    else:
		flag_search=""
	    nintx = npart
	    if psr.p0 > 0.5:
	        nintx = 16
	    elif psr.p0 > 0.1:
	        nintx = 32
	    nbinx = nbin
	    if psr.p0 < 1.7e-3 and nbin > 18:
	        nbinx = 28
	    elif psr.p0 < 1.0e-2 and nbin > 20:
	        nbinx = 50
	    elif psr.p0 < 5.0e-2 and nbin > 30:
	        nbinx = 128
	    prep_str = " -n %d -nsub %d -npart %d -fine%s -mask %s -noxwin %s" \
	        %(nbinx, nsub, nintx, raw_opt, beam_cand.mask[-1], fits)
	    if str(igchn)!="None":
		chans=Popen("grep %s %s | awk '{print $2}'" %(psr.name,igchn),shell=True,stdout=PIPE).communicate()[0].split('\n')[0]
		if len(chans)!=0:
		    prep_str="-ignorechan " +chans+prep_str
	    if nbinx != 200:
	        psr.change_bins(nbinx)		# Mask fraction initially calculated using 50 bins - if this isn't the case, this recalculates with the correct number
	    flag_search = flag_search + psr.specflgs
	    prof_cand = glob('guppi*%s*%s*.pfd.bestprof' %(mjd,beam_cand.num))
	    if len(prof_cand) == 0 and psr.par != '':			# Begin prepfolding, first with a par file (if available), then with command line parameters
	        print "Attempting to fold %s with par file..." %psr.name
	        if len(glob("prepfold_%s_%s_output.txt"%(mjd,beam_cand.num))) > 0:
		    call("rm prepfold_%s_%s_output.txt"%(mjd,beam_cand.num), shell=True)
	        fold_out = open('prepfold_%s_%s_output.txt' %(mjd,beam_cand.num),'w')
		if Popen(["grep PB %s" %psr.par],stdout=PIPE,shell=True).communicate()[0]=='':
	            fold_out.write("prepfold -par %s -nopdsearch %s %s"%(psr.par,flag_search,prep_str))
	            p = Popen("prepfold -par %s -nopdsearch %s %s"%(psr.par,flag_search,prep_str), 
		        shell=True, stdout=fold_out)
		else:
		    fold_out.write("prepfold -par %s %s %s"%(psr.par,flag_search,prep_str))
		    p = Popen("prepfold -par %s %s %s"%(psr.par,flag_search,prep_str),
			shell=True, stdout=fold_out)
	        p.wait()
	        fold_out.flush()
	        fold_out.close()
	        del p
	        prof_cand = glob("guppi*%s*%s*.bestprof" %(mjd,beam_cand.num))
	    if len(prof_cand) == 0:
	        print "Attempting to fold %s with simple parameters..." %psr.name
	        if os.path.isfile("prepfold_%s_%s_output.txt" %(mjd,beam_cand.num)):
		    call("rm prepfold_%s_%s_output.txt" %(mjd,beam_cand.num), shell=True)
	        fold_out = open('prepfold_%s_%s_output.txt' %(mjd,beam_cand.num),'w')
	        prepfold_parameters="prepfold -p %.11f -pd 0.0 -nopdsearch" %psr.p0
		if dmsearch:
		    flag_search = ""
		else:
		    flag_search = "-nodmsearch"
	        if psr.dm != "*":
		    prepfold_parameters=prepfold_parameters+" -dm %.3s " %psr.dm
		prepfold_parameters = prepfold_parameters + psr.specflgs
	        fold_out.write(prepfold_parameters+" %s %s" 
		    %(flag_search, prep_str))
	        p = Popen(prepfold_parameters+" %s %s" 
		    %(flag_search, prep_str), 
		    shell=True,stdout=fold_out)
	        p.wait()
	        fold_out.flush()
	        fold_out.close()
	        del p
	        prof_cand = glob('guppi*%s*%s*.pfd.bestprof' %(mjd,beam_cand.num))
	    if len(prof_cand) == 0:                             # if the prepfolding fails, skip to the next pointing
	        print "Not fully processed, skipping..."
	    else:
	        if dmsearch:                                    # Record the new DM if searching
		    psr.newdm=Popen(["grep 'Best DM' prepfold_%s_%s_output.txt | awk '{print $6}'" %(mjd,beam_cand.num)],stdout=PIPE,shell=True).communicate()[0][0:-1]
		    newdmerr=(psr.p0/nbinx)*(.35**3/100.)*(1/8.298e-6)			    # Calculate error on DM using bin resolution
	        # Deal with multiple bestprof files
	        best_prof_list = np.array(glob('gu*%s*%s*prof' %(mjd,beam_cand.num)))       # Finding the best_prof to use for calculation of flux and S/N
	        if len(best_prof_list) == 1:
		    prof_file = best_prof_list[0]
	        elif len(best_prof_list) > 1:
		    chisq = np.array([float(Popen(["grep 'Reduced' %s | \
		        awk '{print $5}'" %prof_file], stdout=PIPE, shell=True).communicate()[0])
		        for prof_file in best_prof_list])
		    prof_file = best_prof_list[chisq == chisq.max()][0]
	        elif len(best_prof_list) == 0:
		    sys.exit('No bestprof file found for %s' %psr.name)
	        num, val = np.loadtxt(prof_file, dtype='float', unpack=True)          # Read profile information for flux analysis
	        lim = val < val.mean()
		for i in range(3):
		    lim2 = val < val[lim].mean()+2.5*val[lim].std()
		    lim = val < val[lim2].mean()+2.5*val[lim2].std()
		    if val[lim].std() == val[lim2].std():
		        break
	        wren = len(val)-len(val[lim])+0.0
		val=val-val[lim].mean()
		hy=val.max()*0.5
		rts=[]
		if val[0]>hy or val[-1]>hy:
		    val=val.tolist()
		    lim=lim.tolist()
		    lim2=lim2.tolist()
                    val=shift(val,sum(val>hy)+1)
                    lim=shift(lim,sum(val>hy)+1)
		    lim2=shift(lim2,sum(val>hy)+1)
		    val=np.array(val)
		    lim=np.array(lim)
		    lim2=np.array(lim2)
		for i in range(len(val)-1):
		    if i>0 and val[i]>hy:
		        if val[i-1]<hy and val[i+1]<hy:
	                    m=val[i]-val[i-1]
        	            b=val[i]-m*num[i]-hy
                	    rts.append(np.roots([m,b])[0])
                	    m=val[i+1]-val[i]
			    b=val[i]-m*num[i]-hy
			    rts.append(np.roots([m,b])[0])
		        elif val[i-1]<hy:
			    m=val[i]-val[i-1]
			    b=val[i]-m*num[i]-hy
			    rts.append(np.roots([m,b])[0])
		        elif val[i+1]<hy:
			    m=val[i+1]-val[i]
			    b=val[i]-m*num[i]-hy
			    rts.append(np.roots([m,b])[0])
		if len(rts)<2:
		    W=1.0
		elif len(rts)>len(val)*0.75:
		    W=wren
		else:
		    W=0.0
		    for i in range(len(rts)):
		        if i%2==1:
			    W+=rts[i]-rts[i-1]
	        if W == 0.0:							      # If pulse width isn't found, set width to 1 bin (may need to adjust this)
		    W = 1.0
#	        elif W > 1:
#		    high_val = num[val == val.max()][0]
#		    for i in range(len(lim)):					      # If there are single-bin gaps in the pulse width, fill in gaps
#		        if i != high_val and lim[i] == False and \
#			    lim[i-1] == True and lim[(i+1)%len(lim)] == True:
#			    lim[i] = True
#		    W = len(val)-len(val[lim])+0.0
	        P = len(val)+0.0
		#plt.plot(num,val)
		#plt.plot(num[lim],val[lim])
		#for i in rts:
		#    plt.axvline(x=i)
		#plt.show()
	        sig = val[lim].std()
	        if len(beam_cand.fits) != len(beam_cand.maskfrac):					# if Will hasn't reported a mask fraction, calculate it using prepfold output
		    frac_num = float(Popen(["grep points prepfold_%s_%s_output.txt | awk -F'= ' '{print $2}'| head -2 | tail -1" \
		        %(mjd,beam_cand.num)],stdout=PIPE,shell=True).communicate()[0][0:-1])
		    frac_denom =  float(Popen(["grep points prepfold_%s_%s_output.txt | awk -F'= ' '{print $2}'| head -1" \
		        %(mjd,beam_cand.num)],stdout=PIPE,shell=True).communicate()[0][0:-1])
		    beam_cand.maskfrac.append(float(1-frac_num/frac_denom))
		    bandwidth=100.*(1-beam_cand.maskfrac[-1])
		poln=2
	        if psr.s350 == '*':				# Calculate the S/N and flux density for the various pulse widths
		    snr_exp = '*'
	        else:
		    snr_exp = psr.s350*2*np.sqrt(poln*120*bandwidth)*\
		        np.sqrt((P-W)/W)/(psr.temp*np.exp((beam_cand.off/0.3)**2 / 1.5))		
	        if psr.w50b != '*' and psr.s350 != '*':
		    if P>psr.w50b:
		        snr_catW = psr.s350*2*np.sqrt(poln*120*bandwidth)*np.sqrt((P-psr.w50b)/psr.w50b)/(psr.temp*np.exp((beam_cand.off/0.3)**2 / 1.5))
		    else:
		        snr_catW = 'W>P'
	        else: 
		        snr_catW = '*'  
	        if psr.wsmearb != '*' and psr.s350 != '*':
		    if P>psr.wsmearb:
		        snr_smearW = psr.s350*2*np.sqrt(poln*120*bandwidth)*np.sqrt((P-psr.wsmearb)/psr.wsmearb)/(psr.temp*np.exp((beam_cand.off/0.3)**2 / 1.5)) 
		    else:
		        snr_smearW = 'W>P'
	        else: 
		    snr_smearW = '*'  
	        snr_beam = (np.array([(n-val[lim].mean()) for n in
		    val]).sum())/(sig*np.sqrt(W))
	        S_beam = snr_beam*psr.temp*np.sqrt(W/(P-W))/(2*np.sqrt(poln*120*bandwidth))
	        S_offset = np.exp((beam_cand.off/0.3)**2 / 1.5)*S_beam
	        # estimate error in flux using dW=1, dSNR=25%, dT_sys=10K,
	        #  dT_obs=5s, dBW=10MHz
	        # includes estimate of error in angular separation (offset
	        #  from center of beam)
	        dS_offset = np.sqrt((S_beam*np.sqrt((1/4.)**2 + (10/psr.temp)**2 +
		    0.5*(1/24.)**2 + 0.5*(10./bandwidth)**2 +
		    0.5*(P/(W*(P-W)))**2))**2 +
		    (0.1*beam_cand.off*np.exp((beam_cand.off/0.3)**2/1.5))**2)
	        if snr_beam < 0.001:
		    S_offset = 0.0
		    dS_offset = 0.0
	        elif snr_beam < snr_min and snr_beam >= 0.001 and beam_cand.off >= 0.79:
		    S_offset = 100.0
		    dS_offset = 100.0
		out_str="%-12s\t %-5.5s \t %-5.5s \t %-8.8s \t %-.5s " %(psr.name,psr.ra,psr.dec,psr.p0,psr.dm)
	        if dmsearch:								# write output 
		    out_str=out_str+"\t %-.5s(%.2f) " %(psr.newdm,newdmerr)
		out_str=out_str+"\t %-.5s \t %-.5s \t %-.4s \t %-.4s \t %-.5s \t %-.7s " %(beam_cand.off,mjd,bandwidth,psr.spindx,psr.temp,W*psr.p0*1000./P)
		if psr.w50flag:
		    out_str=out_str+"\t %-.4s! \t %-.4s!" %(psr.w50b,psr.wsmearb)
		else:
		    out_str=out_str+"\t %-.4s \t %-.4s" %(psr.w50b,psr.wsmearb)
		if snr_beam < 0.1:
		    snr_beam = 0.0
		out_str=out_str+" \t %-.6s\t %-.6s\t%-.6s\t%-.6s\t   %-.6s    \t %-.6s    \t" %(snr_exp,snr_catW,snr_smearW,snr_beam,S_offset,dS_offset)
		summary_file.write(out_str)
	        if snr_beam > snr_min:				# output to stdout
		    print "PSR %s detected in beam #%s on MJD %d with S/N of %.3f; expected S/N %.6s" \
		        %(psr.name, beam_cand.num, mjd, snr_beam, snr_exp)
		    summary_file.write("Detected in Beam %6s\n" %beam_cand.num)
		    num_detect+=1
		    if psr.name not in name_detect:
		        name_detect.append(psr.name)
		    if not os.path.isdir(work_dir+'detection_plots'):     # Link relevant files in detection plot repository
		        os.chdir(work_dir)
		        call('mkdir detection_plots',shell=True)
		        os.chdir(psr.name+"_temp/")
		    pfd_str = "detection_plots/guppi_%d_GBNCC%s_%s.pfd" \
		        %(mjd,beam_cand.num,psr.name)
		    if len(glob('../%s.*' %pfd_str))==0:
		        if len(glob('g*%s*%s*2bit*ps' %(mjd,beam_cand.num))) > 0:
			    file_str = "%s*%s*2bit" %(mjd,beam_cand.num)
		        else:
			    file_str = "%s_%s" %(mjd,beam_cand.name)
			best_prof_list = glob('g*%s*prof'%file_str)
		        if len(best_prof_list) > 1:
			    check_f = ''
			    check_n = 0.0
			    for prof in best_prof_list:
			        try:
				    cs = float(Popen(["grep 'Reduced' %s | \
				        awk '{print $5}'" %prof_file], stdout=PIPE, shell=True).communicate()[0])
				    if cs > check_n:
				        check_f = prof 
			        except:
				    pass
			    if check_f != '':
			        file_str = check_f[3:-10]
		        os.chdir('%sdetection_plots' %work_dir)
		        call('cp -fs %s%s*/g*%s*ps %s.ps'%(work_dir,psr.name,
			    file_str,pfd_str.split('/')[1]), shell=True)
		        call('cp -fs %s%s*/g*%s*prof %s.bestprof' %(work_dir,
			    psr.name,file_str,pfd_str.split('/')[1]),shell=True)
		        os.chdir("%s%s_temp" %(work_dir,psr.name))
		    if use_comp == 'zuul' and len(glob("%s%s*" %(renee_dir,pfd_str)))==0:
		        os.chdir("%sdetection_plots/" %work_dir)
		        call("ln -fs %s%s_temp/*%s*pfd ." %(work_dir,psr.name,beam_cand.num), shell=True)        
		        os.chdir("%s%s_temp" %(work_dir,psr.name))
		    elif use_comp == 'GB' and len(glob("%s%s*" %(renee_dir,pfd_str))) == 0:
		        call("cp  gu*%s_%s*ps %sdetection_plots/"
			    %(beam_cand.num,psr.name,work_dir), shell=True)
		        call("cp  gu*%s_%s*prof %sdetection_plots/"
			    %(beam_cand.num,psr.name,work_dir), shell=True)
	        else:   					# In the case of a non-detection
		    if snr_exp == '*':
		        exptd = "Unknown exp."
		    elif snr_exp > snr_min*1.5:
		        exptd = "Unexpected - expected S/N %.3f" %snr_exp
		    else:
		        exptd = "Expected"
		    print "PSR %s not detected in beam #%s; S/N < %.1f: %s" \
		        %(psr.name, beam_cand.num, snr_min, exptd)
		    if dmsearch and psr.dm != '*' and abs(float(psr.dm)-float(psr.newdm))>10.0:				# If DM searching and there is a non-detection, print difference in DMs to note drastic changes
		        print "Note: DM changed significantly (difference = %s p/cc)" %(abs(float(psr.dm)-float(psr.newdm)))
		    summary_file.write("Not Detected in %6s\n" %beam_cand.num)
	        os.chdir(work_dir)

print ""        
print "YYYEEEEEAAAAAAA (code ran successfully)"
print ""
name_detect = np.array(name_detect)
if yflg:
    y=[1]

## Summarize results ##
print "#####  Summary  #####"
print "%d pulsars in survey area\nDetected %d pulsars in %d beams\nUnable to find \
observations on %s for %d beams" \
    %(sum(y),len(name_detect),num_detect,use_comp,num_unobs)
summary_file.write("#@\n")
summary_file.write("#####  Summary  #####@\n")
summary_file.write("#@ %d pulsars in survey area\n#@ Detected %d pulsars in %d beams\n#@ Unable to find observations on %s for %d beams" %(sum(y),len(name_detect),num_detect,use_comp,num_unobs))
summary_file.close()
