# (Not-so-)Simple script to find GBNCC beams on zuul and, optionally, process with PRESTO     ###
# Processing requires environment set-up before script execution
# Can make plot of beams (observed, detected, and unobserved) as beamcheck_PSR.png
# Write out information about pulsars found as GBNCC_ATNF_beamcheck_[comp]_[n].txt (or [outrt]_[comp]_[n].txt)
# Combine results from multiple computer systems with combine_beamcheck.py
# Make individual summaries with GBNCC_individual_summaries.py
# Make nice plots of profiles with GBNCC_profile_plots.py
# Written by Renee Spiewak
# Last edit Jun. 5, 2017

# -h,--help   Print help page with flags listed
# -pos        ra and dec in hh:mm:ss.ss (-)dd:mm:ss.ss form
# -p          period in s
# -pd         period derivative in s/s (default: 0)
# -dm         DM in pc/cc
# -f          file name (assumes file has columns of NUMBER, JNAME, RA (deg), DEC (deg), P0 (s), P1 (s/s), DM (cm^-3); include * in unknown ranges)
# -beam       beam number as integer
# -n          number of profile bins (default: 50)
# -nsub       number of channels (default: 128) 
# -npart      number of subintegrations (default: 40)
# -angle      max angular offset between beam and pulsar in deg (default: 0.5 deg)
# -snr        minimum S/N for detection (default: 6) 
# -proc       process files (default: no processing)
# -cen        process central (original) beam
# -comp       select computer system (only 'zuul' or 'GB' allowed)
# -pub        Only process published pulsars
# -rfi        Use Scott's options to remove RFI
# -nonew      Do not copy new fits files into directories 
# -out	      Choose output location
# --par       Provide par file or directory in which par file is stored (if the latter, par file is assumed to be [pulsarJname].par 
# --rew	      Reprocess, i.e. delete previous output
# --outrt     Set root for output file (default is GBNCC_ATNF_beamcheck)
# --dmsearch  Search DMs with prepfold

#########################################################

import numpy as np
import sys, os
from glob import glob
import subprocess as sproc
from time import strftime
import argparse as ap
from astropy.coordinates import SkyCoord,Latitude,Longitude 
from astropy.table import Table,Column
from astropy import units as u, constants as c

# Function to find distance 
def ang_offset(lon1,lat1,lon2,lat2):
	x="%f %f" %(lon1,lat1)
	y="%f %f" %(lon2,lat2)
	return float(SkyCoord(x,unit=("deg","deg")).separation(SkyCoord(y,unit=("deg","deg")))/u.deg)


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
                      help="ra and dec ('hh:mm:ss.ss (-)dd:mm:ss.ss')")
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
    pars.add_argument("-n","--nbin",action="store",default=50,type=int,
                      help="number of profile bins")
    pars.add_argument("--nsub",action="store",default=128,type=int,
                      help="number of channels")
    pars.add_argument("--npart",action="store",default=40,type=int,
                      help="number of subintegrations")
    pars.add_argument("--angle",action="store",default=0.5,type=float,
                      help="max offset between beam and pulsar in deg")
    pars.add_argument("--snr",action="store",default=6,type=float,
                      help="minimum S/N for detection")
    pars.add_argument("--proc",action="store_true",
                      help="process files")
    pars.add_argument("--center",action="store_true",dest="find_psr",
                      help="process central (original) beam")
    pars.add_argument("--pub",action="store_true",
                      help="only process published pulsars")
    pars.add_argument("--rfi",action="store_true",
                      help="use Scott's options to remove RFI")
    pars.add_argument("--nonew",action="store_false",dest="pull_new",
                      help="do not copy new fits files into directories")
    pars.add_argument("-o","--out",action="store",help="output location (default is /users/rspiewak/pulsars/)")
    pars.add_argument("--par",action="store",help="parameter file or directory in which parameter file is stored")
    pars.add_argument("--rew",action="store_true",dest="rewrite",help="delete old output and rewrite")
    pars.add_argument("--outrt",action="store",dest="outputroot",help="root for name of summary file")
    pars.add_argument("--dmsearch",action="store_true",dest="dmsearch",help="search DMs with prepfold")
    args = vars(pars.parse_args())

    if args['pos'] != None and len(args["pos"].split()) > 1:
        ra_astro = Longitude(args["pos"].split()[0],unit="hourangle")
        dec_astro = Latitude(args["pos"].split()[1],unit="deg")
        args["ra_psr"] = ra_astro.value
        args["dec_psr"] = dec_astro.value

    return args

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
		self.beams=[]
		self.par=''
		self.params=name,ra,dec,p0,p1,dm	
		self.pos="%s %s" %(ra,dec)
		if cat_dat[cat_dat['name']==name]['temp'] != '*':
			self.temp=23+float(cat_dat[cat_dat['name']==name]['temp'])
		else:
			self.temp=cat_dat[cat_dat['name']==name]['temp']
		self.w50=cat_dat[cat_dat['name']==name]['W50']
		if self.w50 != '*':
			self.w50b=float(self.w50)/(20.*self.p0)
		else:
			self.w50b='*'
		if dm != '*' and self.w50 != '*':
			tdm=8.3*10**6*float(dm)*df/cent_freq**3
			ts=10**(-6.46+0.154*np.log10(float(dm))+1.07*(np.log10(float(dm)))**2-3.86*np.log10(float(cent_freq)/1000.))
			self.wsmear=np.sqrt(float(self.w50)**2+ts**2+tdm**2)
			self.wsmearb=float(self.wsmear)/(20.*self.p0)
		else:
			self.wsmear='*'
			self.wsmearb='*'
		self.s400=cat_dat[cat_dat['name']==name]['S400'][0]
		self.s1400=cat_dat[cat_dat['name']==name]['S1400'][0]
		self.spindx=cat_dat[cat_dat['name']==name]['SPINDX'][0]
		if self.dec > -40 - ang_max:
			self.north=True
		else:
			self.north=False 
                if self.spindx != '*':
                    if self.s400 != '*':
                        self.s350 = float(self.s400)*(0.35/0.4)**float(self.spindx)
                    elif self.s1400 != '*':
                        self.s350 = float(self.s1400)*(0.35/1.4)**float(self.spindx)
                elif self.s400 != '*' and self.s1400 != '*':
                    self.s350=float(self.s400)**(0.35/0.4)**(np.log(float(self.s400)/float(self.s1400))/-1.253)
                elif self.s400 != '*':
                    self.s350=float(self.s400)*((0.35/0.4)**(-1.7))
                elif self.s1400 != '*':
                    self.s350=float(self.s1400)*((0.35/1.4)**(-1.7))
                else:
                    self.s350='*'


	def add_beam(self,beam):
		self.beams.append(beam)	
	def add_par(self,par_str):
		self.par=par_str
	def add_temp(self,temp):
		self.temp=temp
	def beam_names(self):
		Str=''
                for i in self.beams:   
                	Str=Str+" %s" %i.beams
		print Str
	def beam_numbers(self):
		Str=''
		for i in self.beams:
			Str=Str+" %s" %i.num
		print Str
	def set_pos(self,(ra,dec)):
		self.ra=ra
		self.dec=dec

class Beam:
	def __init__(self,beam):
		self.name=beam
		self.num=beam.strip('GBNCC')
		self.mask=''
                self.maskfrac=''
		self.fits=''
		self.ra=pointings[pointings['pointing']==beam]['RAdeg'][0]
		self.dec=pointings[pointings['pointing']==beam]['Decdeg'][0]
		self.pos="%s %s" %(self.ra,self.dec)
	def add_mask(self,mask):
		self.mask=mask	
                if self.name in mask_fractions['pointing']:
                    self.maskfrac=mask_fractions[mask_fractions['pointing']==self.name]['mask_frac'][0]
        def mask_frac(self,frac):
                self.maskfrac=float(frac)
	def add_fits(self,fits):
		self.fits=fits
		self.mjd=int(self.fits.split('_')[1])
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
pull_new = args["pull_new"]
par_file = args["par"]
rewrite = args["rewrite"]
outrt = args["outputroot"]
if "ra_psr" in args:
    ra_psr = args["ra_psr"]
    dec_psr = args["dec_psr"]
    psr=str(SkyCoord(ra_psr,dec_psr,unit=('hourangle','deg')).to_string('hmsdms').strip('s'))
    psr_name = "J%s%s%s%s" %(psr.split()[0].split('h')[0], psr.split()[0].split('h')[1].split('m')[0],\
	 psr.split()[1].split('d')[0], psr.split()[1].split('d')[1].split('m')[0])


renee_dir = "/users/rspiewak/pulsars/"
cent_freq=350. #MHz
df=0.0244140625  #MHz

if "out" in args:
    work_dir = args["out"]
    if work_dir[-1]!='/':
	work_dir=work_dir+'/'
else:
    work_dir = renee_dir
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
cat_dat=Table.read('/users/amcewen/GBNCC_work/psrcat_v158_full.fits')

num_north=0
beam_rejects = np.array([20051,19907,10451,19906,20312,13691,17237,19064,72172,80142,83626,114632,115242,120582])
pulsar_list=[]

# Check file information (assumes order of data in file)
if not file_psr is None:
    p_coord=[]
    if file_psr.split('.')[-1]!='fits':
    	psr_tbl=Table.read(file_psr,format='ascii')
    else:
	psr_tbl=Table.read(file_psr)
    for i in range(len(psr_tbl)):
	if psr_tbl[i][4]!='*':
	    pulsar_list.append(Pulsar((psr_tbl[i][1],psr_tbl[i][2],psr_tbl[i][3],psr_tbl[i][4],psr_tbl[i][5],psr_tbl[i][6])))
	    num_north+=pulsar_list[-1].north
	else:
	    print "%s omitted for invalid period" %psr_tbl[i][1]
    print "%d/%d pulsars in survey area. Finding closest beams..." %(num_north,len(pulsar_list))
    if num_north==0:
	exit()
    for psr in pulsar_list:
	p_coord.append(SkyCoord(psr.ra,psr.dec,unit=('deg','deg')))
    for i in range(len(pulsar_list)):
	for j in range(len(pointings[p_coord[i].separation(SkyCoord(pointings['RAdeg'],pointings['Decdeg']))<ang_max*u.deg]['pointing'])):
	    if pulsar_list[i].north and not int(pointings[p_coord[i].separation(SkyCoord(pointings['RAdeg'],pointings['Decdeg']))<ang_max*u.deg]['pointing'][j].strip('GBNCC')) in beam_rejects:
		pulsar_list[i].add_beam(Beam(pointings[p_coord[i].separation(SkyCoord(pointings['RAdeg'],pointings['Decdeg']))<ang_max*u.deg]['pointing'][j]))
# If no file is given, read pulsar from input parameters
else:
    pulsar=Pulsar((psr_name,ra_psr,dec_psr,p0_psr,pd_psr,dm_psr))
    num_north+=int(pulsar.north)
    if num_north==1:
	print "%s in survey area. Finding closest beams..." %psr_name
    else:
	print "%s not in survey area." %psr_name
	exit()
    p_coord=SkyCoord(pulsar.ra,pulsar.dec,unit=('hourangle','deg'))
    if str(beam_psr) == "None":
	for j in range(len(pointings[p_coord.separation(SkyCoord(pointings['RAdeg'],pointings['Decdeg']))<ang_max*u.deg]['pointing'])):
	    if pulsar.north and not int(pointings[p_coord.separation(SkyCoord(pointings['RAdeg'],pointings['Decdeg']))<ang_max*u.deg]['pointing'][j].strip('GBNCC')) in beam_rejects:
		pulsar.add_beam(Beam(pointings[p_coord.separation(SkyCoord(pointings['RAdeg'],pointings['Decdeg']))<ang_max*u.deg]['pointing'][j]))
    else:
	if pulsar.north and not int(beam_psr) in beam_rejects:
	    beam_psr='GBNCC'+beam_psr
	    pulsar.add_beam(Beam(beam_psr)) 
    pulsar_list.append(pulsar)

# Create cumulative output file
if "outputroot" in args:
    summary_file = open("%s%s_%s_%dpulsars.txt" %(work_dir,outrt,use_comp,len(pulsar_list)), 'w')
else:
    summary_file = open("%sGBNCC_ATNF_beamcheck_%s_%dpulsars.txt" %(work_dir,use_comp,len(pulsar_list)), 'w')
summary_file.write("# Using S/N > %s for detection threshold on %s\n" %(snr_min,use_comp))
summary_file.write("# Using beams within %.1f deg of pulsars from %s\n" %(ang_max,file_psr))
if '-pub' in args:
    summary_file.write('# Only searched for published pulsars\n')
summary_file.write("# Flux adjusted for offset from beam center\n")
summary_file.write("# S/Ne = Expected S/N estimated from published flux values using width calculated from data; S/Nc = expected S/N calculated using catalog W50; S/Ns = S/N calulated using Weff (smeared out); S/Nm = measured S/N \n")
summary_file.write("# Flux estimated from measured S/N, uncertainty estimated from dependencies\n")
summary_file.write("# File generated %s (ET)\n" %strftime("%Y-%m-%d %H:%M:%S"))
summary_file.write("#   PSRJ  \t RA \t Dec \t P0 \t DM \t Dist  \t MJD     BW\t Pol \t Spindx  Temp \tWdat  \tWcat \tWsmr \t S/Ne \t S/Nc \tS/Ns\tS/Nm \t   Flux\t\tdFlux\t\tStatus \n")
summary_file.write("#\t\t (deg)\t (deg) \t (s) \t (p/cc)\t (deg)\t (day)\t(MHz)\t (num)\t\t (K) \t(bins) \t(bins) \t(bins) \t\t\t\t\t   (mJy)\t(mJy)  \n")

name_detect = []
num_detect = 0
num_unobs = 0

print "Checking S/N for profiles -- threshold = %.1f" %snr_min


# Find observations and, optionally, process with PRESTO
for psr in pulsar_list:
    print ""
    print ""
    print psr.name, psr.beam_numbers()
    print ""
    if psr.name[0] == "C": #or stat_psr == 'L' or stat_psr == 'G':
	unpub_psr = True
    else:
	unpub_psr = False
    if f_pub and unpub_psr:
	continue     # Skip unpublished pulsars when flag used
    if psr.temp == '*':
	psr.add_temp(23+cat_dat[cat_dat['temp']!='*'][np.argsort(SkyCoord(psr.pos,unit=('deg','deg')).separation(SkyCoord(cat_dat[cat_dat['temp']!='*']['RA'],cat_dat[cat_dat['temp']!='*']['Dec'])))[:3]]['temp'].astype('float').mean())
    if psr.dm <= 0 or psr.p0 <= 0:   # Check the DM and P for the current pulsar
	print "PSR %s rejected due to invalid parameters" %psr.name
	continue
    if str(par_file) == "None":
	if len(glob('%sparfiles/%s_short.par' %(renee_dir,psr.name)))>0:
	    psr.add_par(glob('%sparfiles/%s_short.par' %(renee_dir,psr.name))[0])
	if len(glob(data_dir+'amcewen/pars/%s*' %psr.name))>0:
	    psr.add_par(glob(data_dir+'amcewen/pars/%s*' %psr.name)[0])
	if psr.par == '' and len(psr.name)>=10:
	    print "Cannot find par file in %sparfiles/ or in %samcewen/pars/" %(renee_dir,data_dir)
    else:
	if os.path.isdir(par_file):
    	    try:
		psr.add_par(glob(par_file+"/*%s*par" %psr.name)[0])
	    except:
		print "Can't find par file for %s in %s" %(psr.name,par_file)
	else:
	    if os.path.isfile(par_file):
		psr.add_par(par_file)
	    else:
		print "Invalid par file"
    if not os.path.isdir("%s_temp" %psr.name) and pull_new == True:
	os.mkdir("%s_temp" %psr.name)
    elif not os.path.isdir("%s_temp" %psr.name):
	print "Skipping new pulsar: "+psr.name
	continue
    elif rewrite:
	sproc.call("rm -r %s%s_temp" %(work_dir,psr.name),shell=True)
	os.mkdir("%s%s_temp" %(work_dir,psr.name))
    os.chdir("%s_temp/" %psr.name)
    if not os.path.isfile('%s_short.par' %psr.name) and not psr.par == '':
	sproc.call('cp -sf %s .' %psr.par,shell=True)
    for beam_cand in psr.beams:
		print beam_cand.num
		beam_cand.ang_off(psr.ra,psr.dec)
		fits_list = glob('%s*/gu*G????%s_*fits' %(data_dir,beam_cand.num))
		if len(fits_list) > 0:
		    # Process files to find detections
		    if os.getcwd() != work_dir+'%s_temp' %psr.name:
		        if not os.path.isdir(work_dir+'%s_temp' %psr.name): 
			    os.chdir(work_dir)
			    os.mkdir(psr.name+'_temp')
		        os.chdir(work_dir+'%s_temp/' %psr.name)
		    fits_list = glob('%s2*/gu*G????%s_*fits' %(data_dir,beam_cand.num))
		    if len(glob('gu*G????%s_*fits' %beam_cand.num))==0 and len(fits_list)>0 and pull_new==True:
		        sproc.call("cp -s %s ." %fits_list[-1], shell=True)
		    elif len(glob('gu*G????%s_*fits' %beam_cand.num))==0 and len(fits_list)>0:
		        print "Skipping new data"
		        continue
		    A = sproc.Popen(["ls gu*G????%s_*fits" %beam_cand.num], stdout=sproc.PIPE, shell=True)
		    beam_cand.add_fits(A.communicate()[0][0:-1])
		    rfi_std = "-time 1 -timesig 2 -freqsig 3"
		    if "2bit" in beam_cand.fits and rfi_fil == True:
		        raw_opt = " -noscales -noweights"
		        rfi_std = "-time 1 -freqsig 3"
		    else:
		        raw_opt = " "
		    if beam_cand.mjd > 56710 and rfi_fil == True:
		        rfi_opt = " -zapchan 2470:3270"
		    else:
		        rfi_opt = " "
		    if len(glob('%s%s_temp/*%s*%s*rfifind.mask' %(work_dir,psr.name,beam_cand.mjd,beam_cand.num)))==0 \
			    or len(glob('%s%s_temp/*%s*%s*rfifind.stats' %(work_dir,psr.name,beam_cand.mjd,beam_cand.num)))==0:
		        if len(glob(data_dir+'amcewen/mask_files/*%s*%s*rfifind.mask' %(beam_cand.mjd,beam_cand.num)))==0 \
			    or len(glob(data_dir+'amcewen/mask_files/*%s*%s*rfifind.stats' %(beam_cand.mjd,beam_cand.num)))==0:
			        mask_dir=glob(data_dir+'20*/*%s*GBNCC%s*' %(beam_cand.mjd,beam_cand.num))[0].split('/')[5] 
			        if len(glob(data_dir+mask_dir+'/*rfi*tar*'))>0:
				    print "untarring mask file"
				    os.chdir(data_dir+'amcewen/mask_files/')
				    for tar_file in glob(data_dir+mask_dir+'/*rfi*tar*'):
					    if len(glob('*%s*%s*rfifind.mask' %(beam_cand.mjd,beam_cand.num)))==0:
						    sproc.call('tar -xf '+tar_file+' \*GBNCC\*d.mask', shell=True)
					    if len(glob('*%s*%s*rfifind.stats' %(beam_cand.mjd,beam_cand.num)))==0:
						    sproc.call('tar -xf '+tar_file+' \*GBNCC\*d.stats', shell=True)
				    sproc.call('cp *%s*GBNCC%s*rfifind.mask %s*%s*/' %(beam_cand.mjd,beam_cand.num,work_dir,psr.name),shell=True)
				    sproc.call('cp *%s*GBNCC%s*rfifind.stats %s*%s*/' %(beam_cand.mjd,beam_cand.num,work_dir,psr.name),shell=True)
				    os.chdir(work_dir+'%s_temp' %psr.name)
			        else:
				    if len(glob('*_%s_%s_rfifind.stats' %(beam_cand.mjd,beam_cand.num))) == 0\
					    or len(glob('*_%s_%s_rfifind.mask' %(beam_cand.mjd,beam_cand.num))) == 0:
				        print "Running rfifind for %s beam candidate %s..." %(psr.name,beam_cand.num)
				        rfi_out = open('rfifind_%s_%s_output.txt' %(beam_cand.mjd,beam_cand.num),'w')
				        rfi_out.write("nice rfifind %s%s%s -o %s_%s_%s %s"
					          %(rfi_std,raw_opt,rfi_opt,psr.name,beam_cand.mjd,beam_cand.num,beam_cand.fits))
				        p = sproc.Popen("nice rfifind %s%s%s -o %s_%s_%s %s"
						    %(rfi_std,raw_opt,rfi_opt,psr.name,beam_cand.mjd,beam_cand.num,beam_cand.fits),
						    stdout=rfi_out, shell=True)
				        p.wait()
				        rfi_out.flush()
				        rfi_out.close()
				        del p
				        sproc.call('rm *.bytemask *.inf *d.ps *.rfi', shell=True)
			        #print "retrieving .mask and .stats files from %samcewen/mask_files/" %data_dir
			        sproc.call('cp *%s*%s_*d.mask %samcewen/mask_files/' %(beam_cand.mjd,beam_cand.num,data_dir), shell=True)
			        sproc.call('cp *%s*%s_*d.stats %samcewen/mask_files/' %(beam_cand.mjd,beam_cand.num,data_dir), shell=True) 
		        else:
			        #print "retrieving .mask & .stats files from %samcewen/mask_files" %data_dir
			        mask_dir=glob(data_dir+'amcewen/mask_files/*%s*%s_*rfifind.mask' %(beam_cand.mjd,beam_cand.num))[0]
			        stats_dir=glob(data_dir+'amcewen/mask_files/*%s*%s_*rfifind.stats' %(beam_cand.mjd,beam_cand.num))[0]
			        sproc.call('cp -fs %s .' %mask_dir,shell=True)
			        sproc.call('cp -fs %s .' %stats_dir,shell=True)
		    else:
		        #print "retrieving .mask & .stats files from %s%s_temp/" %(work_dir,psr.name)
		        mask_list=glob(work_dir+'%s_temp/*%s*%s*rfifind.mask' %(psr.name,beam_cand.mjd,beam_cand.num))
		        stats_list=glob(work_dir+'%s_temp/*%s*%s*rfifind.stats' %(psr.name,beam_cand.mjd,beam_cand.num))
		        for mask_dir in mask_list:
			    stats_dir=mask_dir.strip('mask')+'stats' 
			    if os.getcwd().split('/')[-1] != mask_dir.split('/')[-2]: 
			        sproc.call('cp -fs %s .' %mask_dir, shell=True)
			        sproc.call('cp -fs %s .' %stats_dir, shell=True)
		    #print "mask file located"
		    mask_file=glob('*%s*%s*rfifind.mask' %(beam_cand.mjd,beam_cand.num))[0]				
		    beam_cand.add_mask(mask_file)
		    if beam_cand.maskfrac != '':
		        bandwidth=100.*(1-beam_cand.maskfrac)
		    if "dmsearch" not in args:
		        flag_search = '-nosearch'
		        if (len(psr.name) < 10 and psr.p0 < 0.1):
			    flag_search = '-nodmsearch'
		    elif len(psr.name) < 10 and psr.p0 < 0.1:
		        flag_search=""
		    else:
		        flag_search="-nopsearch"
		    nintx = 64
		    if psr.p0 > 0.5:
		        nintx = 16
		    elif (psr.p0 > 0.1) or (psr.p0 > 0.5):
		        nintx = 32
		    nbinx = nbin
		    if psr.p0 < 1.7e-3 and nbin > 18:
		        nbinx = 16
		    elif psr.p0 < 2.5e-3 and nbin > 20:
		        nbinx = 20
		    elif psr.p0 < 3.5e-3 and nbin > 30:
		        nbinx = 30
		    prep_str = "-n %d -nsub %d -npart %d -fine%s -mask %s -noxwin %s" \
			       %(nbinx, nsub, nintx, raw_opt, beam_cand.mask, beam_cand.fits)
		    prof_cand = glob('guppi*%s*%s*.pfd.bestprof' %(beam_cand.mjd,beam_cand.num))
		    if len(prof_cand) == 0 and not psr.par == '':
		        print "Attempting to fold %s with par file..." %psr.name
		        if len(glob("prepfold_%s_%s_output.txt"%(beam_cand.mjd,beam_cand.num))) > 0:
			    sproc.call("rm prepfold_%s_%s_output.txt"%(beam_cand.mjd,beam_cand.num), shell=True)
		        fold_out = open('prepfold_%s_%s_output.txt' %(beam_cand.mjd,beam_cand.num),'w')
		        fold_out.write("prepfold -par %s %s %s"%(psr.par,flag_search,prep_str))
		        p = sproc.Popen("prepfold -par %s %s %s"%(psr.par,flag_search,prep_str), 
				        shell=True, stdout=fold_out)
		        p.wait()
		        fold_out.flush()
		        fold_out.close()
		        del p
		        prof_cand = glob("guppi*%s*%s*.bestprof" %(beam_cand.mjd,beam_cand.num))
		    if len(prof_cand) == 0:
		        print "Attempting to fold %s with simple parameters..." %psr.name
		        if os.path.isfile("prepfold_%s_%s_output.txt" %(beam_cand.mjd,beam_cand.num)):
			    sproc.call("rm prepfold_%s_%s_output.txt" %(beam_cand.mjd,beam_cand.num), shell=True)
		        fold_out = open('prepfold_%s_%s_output.txt' %(beam_cand.mjd,beam_cand.num),'w')
		        prepfold_parameters="prepfold -p %.11f " %psr.p0
		        if psr.p1 != "*":
			    prepfold_parameters=prepfold_parameters+"-pd %s " %psr.p1
		        if psr.dm != "*":
			    prepfold_parameters=prepfold_parameters+"-dm %.3s " %psr.dm
		        fold_out.write(prepfold_parameters+" %s %s" 
				       %(flag_search, prep_str))
		        p = sproc.Popen(prepfold_parameters+" %s %s" 
				        %(flag_search, prep_str), 
				        shell=True,stdout=fold_out)
		        p.wait()
		        fold_out.flush()
		        fold_out.close()
		        del p
		        prof_cand = glob('guppi*%s*%s*.pfd.bestprof' %(beam_cand.mjd,beam_cand.num))
		    if len(prof_cand) == 0:
		        print "Not fully processed, skipping..."
		    else:
		        # Cut extra RFI "manually"
		        if int(beam_cand.fits.split("_")[1]) > 56710 and rfi_fil == True \
		           and len(glob('prepfold_%s_*txt')) > 0:
			    p = sproc.Popen("tail -4 prepfold_%s_%s_*txt | head -1" %(beam_cand.mjd,beam_cand.num),
					    shell=True,stdout=sproc.PIPE).communicate()[0]
			    if float(p.split()[-1]) > 1.3:
			        print "Removing extra RFI"
			        pfd_out = open("show_pfd_%s_%s_output.txt" %(beam_cand.mjd,beam_cand.num),'w')
			        p = sproc.Popen('show_pfd -killsubs 63,77:102 -noxwin %s*pfd'
					    %beam_cand.fits.split('.')[0],shell=True,stdout=pfd_out)
			        p.wait()
			        pfd_out.flush()
			        pfd_out.close()
			        del p
		        # Deal with multiple bestprof files
		        best_prof_list = np.array(glob('gu*%s*%s*prof' %(beam_cand.mjd,beam_cand.num))) 
		        if len(best_prof_list) == 1:
			    prof_file = best_prof_list[0]
		        elif len(best_prof_list) > 1:
			    chisq = np.array([float(sproc.Popen(["grep 'Reduced' %s | \
    awk '{print $5}'" %prof_file], stdout=sproc.PIPE, shell=True).communicate()[0])
					      for prof_file in best_prof_list])
			    prof_file = best_prof_list[chisq == chisq.max()][0]
		        elif len(best_prof_list) == 0:
			    sys.exit('No bestprof file found for %s' %psr.name)
		        num, val = np.loadtxt(prof_file, dtype='float', unpack=True)
		        lim = val < val.mean()
		        i = 0
		        while i < 3:
			    lim2 = val < val[lim].mean()+2.5*val[lim].std()
			    lim = val < val[lim2].mean()+2.5*val[lim2].std()
			    if val[lim].std() == val[lim2].std():
			        break
			    i = i+1
		        W = len(val)-len(val[lim])+0.0
		        if W == 0.0:
			    W = 1.0
		        elif W > 1:
			    high_val = num[val == val.max()][0]
			    for i in range(len(lim)):
			        if i != high_val and lim[i] == False and \
			           lim[i-1] == True and lim[(i+1)%len(lim)] == True:
				    lim[i] = True
			    W = len(val)-len(val[lim])+0.0
		    	
		        P = len(val)+0.0
		        sig = val[lim].std()
		        if beam_cand.maskfrac == '':
			    frac_num = float(sproc.Popen(["cat prepfold_%s_%s_output.txt | grep points | awk -F'= ' '{print $2}'| head -2 | tail -1" \
			        %(beam_cand.mjd,beam_cand.num)],stdout=sproc.PIPE,shell=True).communicate()[0][0:-1])
			    frac_denom =  float(sproc.Popen(["cat prepfold_%s_%s_output.txt | grep points | awk -F'= ' '{print $2}'| head -1" \
			        %(beam_cand.mjd,beam_cand.num)],stdout=sproc.PIPE,shell=True).communicate()[0][0:-1])
			    beam_cand.mask_frac(1-frac_num/frac_denom)
			    bandwidth=100.*(1-beam_cand.maskfrac)
		        poln=float(sproc.Popen(["readfile %s | grep poln | awk '{print $5}'" %beam_cand.fits],stdout=sproc.PIPE,shell=True).communicate()[0][0])
		        if psr.s350 == '*':
			    snr_exp = '*'
		        else:
			    snr_exp = psr.s350*2*np.sqrt(poln*120*bandwidth)*\
				      np.sqrt((P-W)/W)/(psr.temp*np.exp((beam_cand.off/0.3)**2 / 1.5))
			if psr.w50b != '*':
				snr_catW = psr.s350*2*np.sqrt(poln*120*bandwidth)*np.sqrt((P-psr.w50b)/psr.w50b)/(psr.temp*np.exp((beam_cand.off/0.3)**2 / 1.5))
			else: 
				snr_catW = '*'  
			if psr.wsmearb != '*':
				try:
				    snr_smearW = psr.s350*2*np.sqrt(poln*120*bandwidth)*np.sqrt((P-psr.wsmearb)/psr.wsmearb)/(psr.temp*np.exp((beam_cand.off/0.3)**2 / 1.5)) 
				except:
				    snr_smearW = '>P'
			else: 
				snr_smearW = '*'  
		        snr_beam = np.abs(np.array([(n-val[lim].mean()) for n in
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
		        summary_file.write("%-12s\t %-5.5s \t %-5.5s \t %-.5s \t %-.5s \t %-.5s \t %-.5s \t %-.4s \t %-.1s \t %-.4s \t %-.5s \t %-.3s \t %-.3s \t %-.3s \t %-.6s\t %-.6s\t%-.6s\t%-.6s\t   %-.6s    \t %-.6s    \t" \
			    %(psr.name,psr.ra,psr.dec,psr.p0,psr.dm,beam_cand.off,beam_cand.mjd,bandwidth,poln,psr.spindx,psr.temp,W,psr.w50b,psr.wsmearb,snr_exp,snr_catW,snr_smearW,snr_beam,S_offset,dS_offset))
		        if snr_beam > snr_min:
			    print "PSR %s detected in beam #%s on MJD %d with S/N of %.3f; expected S/N %.6s" \
			        %(psr.name, beam_cand.num, beam_cand.mjd, snr_beam, snr_exp)
			    summary_file.write("Detected in Beam #%6s\n" %beam_cand.num)
			    num_detect+=1
			    if psr.name not in name_detect:
			        name_detect.append(psr.name)
			    if not os.path.isdir(work_dir+'detection_plots'):
			        os.chdir(work_dir)
			        sproc.call('mkdir detection_plots',shell=True)
			        os.chdir(psr.name+"_temp/")
			    pfd_str = "detection_plots/guppi_%d_GBNCC%s_%s.pfd" \
				      %(beam_cand.mjd,beam_cand.num,psr.name)
			    if len(glob('../%s.*' %pfd_str))==0:
			        if len(glob('g*%s*%s*2bit*ps' %(beam_cand.mjd,beam_cand.num))) > 0:
				    file_str = "%s*%s*2bit" %(beam_cand.mjd,beam_cand.num)
			        else:
				    file_str = beam_cand.num
			        if len(glob('g*%s*ps'%file_str)) > 1:
				    best_prof_list = glob('g*%s*prof'%file_str)
				    check_f = ''
				    check_n = 0.0
				    for prof in best_prof_list:
				        try:
					    cs = float(sproc.Popen(["grep 'Reduced' %s | \
    awk '{print $5}'" %prof_file], stdout=sproc.PIPE, shell=True).communicate()[0])
					    if cs > check_n:
					        check_f = prof 
				        except:
					    pass
				    if check_f != '':
				        file_str = check_f[3:-10]
			        os.chdir('%sdetection_plots' %work_dir)
			        sproc.call('cp -fs %s%s*/g*%s*ps %s.ps'%(work_dir,psr.name,
					    file_str,pfd_str.split('/')[1]), shell=True)
			        sproc.call('cp -fs %s%s*/g*%s*prof %s.bestprof' %(work_dir,
					    psr.name,file_str,pfd_str.split('/')[1]),shell=True)
			        os.chdir("%s%s_temp" %(work_dir,psr.name))
			    if use_comp == 'zuul' and len(glob("%s%s*" %(renee_dir,pfd_str)))==0:
			        os.chdir("%sdetection_plots/" %work_dir)
			        sproc.call("ln -fs %s%s_temp/*%s*pfd ." %(work_dir,psr.name,beam_cand.num), shell=True)        
			        os.chdir("%s%s_temp" %(work_dir,psr.name))
			    elif use_comp == 'GB' and len(glob("%s%s*" %(renee_dir,pfd_str))) == 0:
			        sproc.call("cp  gu*%s_%s*ps %sdetection_plots/"
				           %(beam_cand.num,psr.name,work_dir), shell=True)
			        sproc.call("cp  gu*%s_%s*prof %sdetection_plots/"
				           %(beam_cand.num,psr.name,work_dir), shell=True)
		        else:
			    if snr_exp == '*':
			        exptd = "Unknown exp."
			    elif snr_exp > snr_min*1.5:
			        exptd = "Unexpected - expected S/N %.3f" %snr_exp
			    else:
			        exptd = "Expected"
			    print "PSR %s not detected in beam #%s; S/N < %.1f: %s" \
			        %(psr.name, beam_cand.num, snr_min, exptd)
			    summary_file.write("Not Detected in #%6s\n" %beam_cand.num)
		        os.chdir(work_dir)
		else:
		    print "Cannot find beam #%s on %s for PSR %s" %(beam_cand.num, use_comp, psr.name)
		    
		    summary_file.write("%-12s\t %-5.5s \t %-5.5s \t %-.5s \t %-.5s \t %-.5s \t * \t * \t * \t %-.4s \t %-.5s \t * \t %-.3s \t %-.3s \t * \t * \t* \t * \t   * \t\t * \t\t Beam #%6s Not Found\n" \
			%(psr.name,psr.ra,psr.dec,psr.p0,psr.dm,beam_cand.off,psr.spindx,psr.temp,psr.w50b,psr.wsmearb,beam_cand.num))
		    num_unobs += 1
		    os.chdir(work_dir)




print ""        
print "YYYEEEEEAAAAAAA (code ran successfully)"
print ""
name_detect = np.array(name_detect)


## Summarize results ##
print "#####  Summary  #####"
print "%d pulsars in survey area\nDetected %d pulsars in %d beams\nUnable to find \
observations on %s for %d beams" \
    %(num_north,len(name_detect),num_detect,use_comp,num_unobs)

