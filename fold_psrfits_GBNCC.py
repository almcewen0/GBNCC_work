# -h,--help   Print help page with flags listed
# --pos       ra and dec in hh:mm:ss.ss (-)dd:mm:ss.ss form
# -p          period in s
# -pd         period derivative in s/s (default: 0)
# --snr       minimum S/N for detection (default: 6) 
# --dm        DM in pc/cc
# -f          file name (assumes file has columns of NUMBER, JNAME, RA (deg), DEC (deg), P0 (s), P1 (s/s), DM (cm^-3); include * in unknown ranges)
# --beam      beam number as integer
# --angle     max angular offset between beam and pulsar in deg (default: 0.5 deg)
# --center    process central (original) beam
# --comp      select computer system (only 'zuul' or 'GB' allowed)
# --rfi       Use Scott's options to remove RFI
# -nonew      Do not copy new fits files into directories 
# -o	      Choose output location (default is /lustre/cv/projects/GBNCC/amcewen)
# --par       Provide par file or directory in which par file is stored (if the latter, par file is assumed to be [pulsarJname].par) 
# --rew	      Reprocess, i.e. delete previous output
# --outrt     Set root for output file (default is fold_psrfits_GBNCC)
# --beamfile  Provide file of beams with pulsar names to search only those beams

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
    pars = ap.ArgumentParser(description="Locate and process data on known pulsars")
    pars.add_argument("--comp",action="store",choices=["zuul","GB"],required=True,help="select computer system")
    pars.add_argument("--pos",action="store",help="ra and dec ('hms dms' or degrees)")
    pars.add_argument("-p","--period",action="store",default=1.0,type=float,help="period in s")
    pars.add_argument("-pd","--pdot",action="store",default=0,type=float,help="period derivative in s/s")
    pars.add_argument("--dm",action="store",default=50,type=float,help="DM in pc/cc")
    pars.add_argument("-f","--file",action="store",help="file name (assumes file has columns of NUMBER, JNAME, RA, DEC, P0, P1, DM)")
    pars.add_argument("--beam",action="store",help="beam number as integer")
    pars.add_argument("--snr",action="store",default=6,type=float,help="minimum S/N for detection")
    pars.add_argument("--angle",action="store",default=0.5,type=float,help="max offset between beam and pulsar in deg")
    pars.add_argument("--center",action="store_true",dest="find_psr",help="process only the closest beam")
    pars.add_argument("--rfi",action="store_true",help="use Scott's options to remove RFI")
    pars.add_argument("--nonew",action="store_true",dest="nonew",help="do not search for .fits files")
    pars.add_argument("-o","--out",action="store",help="output location (default is /lustre/cv/projects/GBNCC/amcewen/)")
    pars.add_argument("--par",action="store",help="parameter file or directory in which parameter file is stored")
    pars.add_argument("--rew",action="store_true",dest="rewrite",help="delete old output and rewrite")
    pars.add_argument("--outrt",action="store",dest="outputroot",help="root for name of summary file")
    pars.add_argument("--beamfile",action="store",dest="beamfile",help="file with list of pulsar names and beam candidates to search")
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
		self.beams=[]
		self.par=''
		self.params=name,ra,dec,p0,p1,dm	
		self.pos="%s %s" %(ra,dec)
		if self.dec > -40 - ang_max:
			self.north=True
		else:
			self.north=False 

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
		self.fits=[]
		self.mjd=[]
		self.ra=pointings[pointings['pointing']==beam]['RAdeg'][0]
		self.dec=pointings[pointings['pointing']==beam]['Decdeg'][0]
		self.pos="%s %s" %(self.ra,self.dec)
        def add_fits(self,fits):
                self.fits.append(fits)
                self.mjd.append(int(fits.split('_')[1]))
	def ang_off(self,ra,dec):
		self.off=ang_offset(self.ra,self.dec,ra,dec)

# Begin full program 
args = proc_args()
yflg=False

p0_psr = args["period"]
pd_psr = args["pdot"]
dm_psr = args["dm"]
file_psr = args["file"]
beam_psr = args["beam"]
ang_max = args["angle"]
snr_min = args["snr"]
find_psr = args["find_psr"]
use_comp = args["comp"]
rfi_fil = args["rfi"]
par_file = args["par"]
rewrite = args["rewrite"]
outrt = args["outputroot"]
nonew = args["nonew"]
beamfile=args["beamfile"]

# Create name of pulsar from command line input
if "ra_psr" in args:
    ra_psr = args["ra_psr"]
    dec_psr = args["dec_psr"]
    psr=str(SkyCoord(ra_psr,dec_psr,unit=('hourangle','deg')).to_string('hmsdms').strip('s'))
    psr_name = "J%s%s%s%s" %(psr.split()[0].split('h')[0], psr.split()[0].split('h')[1].split('m')[0],\
	 psr.split()[1].split('d')[0], psr.split()[1].split('d')[1].split('m')[0])

# Survey specs and directory structure
renee_dir = "/users/rspiewak/pulsars/"

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
#if not os.path.isfile('./GBNCC_pointings.fits'):
#	filename='/lustre/cv/projects/GBNCC/amcewen/GBNCC_posns_by_dec_ALLGBTSKY.txt'
#	np.fromfile(filename)
#	t=Table.read(filename,format='ascii',names=('pointing','RA','Dec'))
#	s=SkyCoord(pointings['RA'],pointings['Dec'],unit=('deg','deg'))
#	t.add_column(Column(s.ra,name='RAdeg'))
#	t.add_column(Column(s.dec,name='Decdeg'))
#	t.write('GBNCC_pointings.fits')
pointings=Table.read('/users/amcewen/GBNCC_work/GBNCC_pointings.fits')
mask_fractions=Table.read('/users/amcewen/GBNCC_work/GBNCC_mask_fraction.fits')
cat_dat=Table.read('/users/amcewen/GBNCC_work/psrcat_v159_full.fits')

beam_rejects = np.array([20051,19907,10451,19906,20312,13691,17237,19064,72172,80142,83626,114632,115242,120582])

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
    else:
	print "%s not valid and/or in survey area." %psr_name
	exit()

print "Checking S/N for profiles -- threshold = %.1f" %snr_min

# Find observations and, optionally, process with PRESTO
for i in range(len(psr_tbl[y])):
    if yflg:
	psr=pulsar
    else:
        psr=Pulsar((psr_tbl[y][i][1],psr_tbl[y][i][2],psr_tbl[y][i][3],psr_tbl[y][i][4],psr_tbl[y][i][5],psr_tbl[y][i][6]))
    p_coord=SkyCoord(psr.ra,psr.dec,unit=('deg','deg'))
    if str(beam_psr) == "None":
	if not find_psr and str(beamfile) == "None":
            for j in range(len(pointings[p_coord.separation(SkyCoord(pointings['RAdeg'],pointings['Decdeg']))<ang_max*u.deg]['pointing'])):
                if not int(pointings[p_coord.separation(SkyCoord(pointings['RAdeg'],pointings['Decdeg']))<ang_max*u.deg]['pointing'][j].strip('GBNCC')) in beam_rejects:
                    psr.beams.append(Beam(pointings[p_coord.separation(SkyCoord(pointings['RAdeg'],pointings['Decdeg']))<ang_max*u.deg]['pointing'][j]))
	elif not find_psr:
	    for beam in sproc.Popen(["grep %s %s" %(psr.name,beamfile)],stdout=sproc.PIPE,shell=True).communicate()[0].split()[1:]:
	        psr.beams.append(Beam('GBNCC'+beam))
	else: 
	    psr.beams.append(Beam(pointings[p_coord.separation(SkyCoord(pointings['RAdeg'],pointings['Decdeg']))==p_coord.separation(SkyCoord(pointings['RAdeg'],pointings['Decdeg'])).min()][0][0]))
    elif int(beam_psr) not in beam_rejects:
	beam_psr='GBNCC'+beam_psr
	psr.beams.append(Beam(beam_psr))
    print ""
    print ""
    print psr.name
    psr.beam_numbers()
    if str(par_file) == "None":   # Unless a .par file/directory is given, search through repositories for match
	if len(glob('%sparfiles/%s_short.par' %(renee_dir,psr.name)))>0:
	    psr.par=glob('%sparfiles/%s_short.par' %(renee_dir,psr.name))[0]
	if len(glob(data_dir+'amcewen/pars/%s*' %psr.name))>0:
	    psr.par=glob('/lustre/cv/projects/GBNCC/amcewen/pars/%s*' %psr.name)[0]
	if psr.par == '':
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
    elif rewrite:
	sproc.call("rm -r %s%s_temp" %(work_dir,psr.name),shell=True)
	sproc.call("rm %sdetection_plots/*" %work_dir,shell=True)
	os.mkdir("%s%s_temp" %(work_dir,psr.name))
    os.chdir("%s_temp/" %psr.name)
    if not os.path.isfile('%s_short.par' %psr.name) and not psr.par == '':      # Linking .par file to directory
	sproc.call('cp -sf %s .' %psr.par,shell=True)
    for beam_cand in psr.beams:
	if os.getcwd().split('/')[-1] != "%s_temp" %psr.name:
	    os.chdir(work_dir+"%s_temp/" %psr.name)
	beam_cand.ang_off(psr.ra,psr.dec)              # Calculate the angular offset between pulsar and beam and store this value in the beam metadata
	if len(glob('*%s*fits' %beam_cand.name))==0:
	    if not nonew:
	        fits_list=glob('%s2*/gu*G????%s_*fits' %(data_dir,beam_cand.num))  # Search for GBNCC .fits files
	        if len(fits_list)==1:
	            sproc.call("cp %s ." %fits_list[0],shell=True)
		    beam_cand.add_fits(fits_list[0].split('/')[-1])       
	        elif len(fits_list)>1:
		    print "Multiple files for pointing #%s" %beam_cand.num
		    for i in fits_list:
	                sproc.call("cp %s ." %i,shell=True)
		        beam_cand.add_fits(i.split('/')[-1])
                else:
                    print "Cannot find beam #%s on %s for PSR %s" %(beam_cand.num, use_comp, psr.name)
                    os.chdir(work_dir)
                    continue
	    else:
		print "Cannot find beam #%s on %s for PSR %s" %(beam_cand.num, use_comp, psr.name)
		os.chdir(work_dir)
		continue
	else:
	    for i in glob('guppi*%s*fits' %beam_cand.name):
	        beam_cand.add_fits(i)
	# Process files to find detections
	for fits,mjd in zip(beam_cand.fits,beam_cand.mjd):
	    num=fits.split('_')[2].strip('GBNCC')
	    if os.getcwd().split('/')[-1] != "%s_temp" %psr.name:
		os.chdir("%s_temp/" %psr.name)
	   # sproc.call('fold_psrfits -P %s -t 10 %s' %(psr.par, fits),shell=True)
	    sproc.call('fits_modhead %s BACKEND GUPPI2bit' %fits,shell=True)
	    sproc.call('dspsr -s -scloffs -E %s -A -e ar -b 256 -O SP_%s_%s_GBNCC %s' %(psr.par,psr.name,num,fits),shell=True)
