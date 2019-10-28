import numpy as np
from astropy.table import Table
from astropy.coordinates import SkyCoord
import subprocess as sproc
from glob import glob
import os
import argparse as ap


class Pulsar:
        def __init__(self,(name,ra,dec,p0,p1,dm)):
                self.name=name
                try:
                        ra.split(':')[1]
                        pos=SkyCoord(ra,dec,unit=('hourangle','deg'))
                except:
                        pos=SkyCoord(ra,dec,unit=('deg','deg'))
                self.ra=pos.ra.value
                self.dec=pos.dec.value
                self.p0=float(p0)
                self.p1=p1
                self.dm=dm
                self.newdm=''
                self.beams=[]
                self.par=''
                self.params=name,ra,dec,p0,p1,dm
                self.pos="%s %s" %(ra,dec)
                self.temp=tskypy(pos.transform_to('galactic').b.value,pos.transform_to('galactic').l.value)
        def beam_numbers(self):
                Str=''
                for i in self.beams:
                        Str=Str+"%s " %i.num
                print Str

class Beam:
        def __init__(self,beam):
                self.name=beam
		try:
		    int(beam)
		    self.num=beam
		    self.name=beam_prefix+str(beam)
		except:
		    self.num=beam.strip(beam_prefix)
                self.mask=[]
                self.maskfrac=[]
                self.fits=[]
                self.mjd=[]
                self.ra=pointings[pointings['beam']==beam]['rajd'][0]
                self.dec=pointings[pointings['beam']==beam]['decjd'][0]
                self.pos="%s %s" %(self.ra,self.dec)
        def add_fits(self,fits):
                self.fits.append(fits)
                self.mjd.append(int(fits.split('_')[1]))
        def ang_off(self,ra,dec):
                self.off=ang_offset(self.ra,self.dec,ra,dec)

def dot(l1,l2):
      l1=l1.astype(int)
      l2=l2.astype(int)
      if len(l1)==len(l2):
          l=l1+l2-1
          for i in range(len(l)):
             if l[i]!=1:
                l[i]=0
      return l.astype(bool)

def tskypy(gb,gl):
    """ Calculate tsky from Haslam table, scale to survey frequency"""
    tskypath = "/users/amcewen/GBNCC_work/tsky.ascii"
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
    # scale temperature and add other sources of heat (t_notsky) before returning
    t_notsky=23
    return tsky_haslam * (cent_freq/408.0)**(-2.6)+t_notsky

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
    pars.add_argument("-p","--period",action="store",default=1.0,type=float,                    help="period in s")
    pars.add_argument("-pd","--pdot",action="store",default=0,type=float,                       help="period derivative in s/s")
    pars.add_argument("--dm",action="store",default=50,type=float,                              help="DM in pc/cc")
    pars.add_argument("-f","--file",action="store",                                             help="file name (assumes file has columns of NUMBER, JNAME, RA, DEC, P0, P1, DM)")
    pars.add_argument("--beam",action="store",                                                  help="beam number as integer")
    pars.add_argument("-n","--nbin",action="store",type=int,                        		help="number of profile bins")
    pars.add_argument("--nsub",action="store",default=128,type=int,                             help="number of frequency channels")
    pars.add_argument("--npart",action="store",default=64,type=int,                             help="number of subintegrations")
    pars.add_argument("--angle",action="store",default=0.5,type=float,                          help="max offset between beam and pulsar in deg")
    pars.add_argument("--snr",action="store",default=6,type=float,                              help="minimum S/N for detection")
    pars.add_argument("--closest",action="store_true",	                          		help="process only the closest beam")
    pars.add_argument("--dat",action="store",default='/lustre/cv/projects/GBNCC/',              help="directory where folders of data are located")
    pars.add_argument("--point",action="store",required=True,                                   help="file with survey pointings (assumes beam,rajd,decjd)")
    pars.add_argument("--rfi",action="store_true",                                              help="use Scott's options to remove RFI")
    pars.add_argument("--par",action="store",                                                   help="parameter file or directory in which parameter file is stored")
    pars.add_argument("--rew",action="store_true",dest="rewrite",                               help="delete old output and rewrite")
    pars.add_argument("--outrt",action="store",dest="outputroot",                               help="root for name of summary file")
    pars.add_argument("-o","--out",action="store",default='/lustre/cv/projects/GBNCC/amcewen/'  help="output location (default is /lustre/cv/projects/GBNCC/amcewen/)")
    pars.add_argument("--dmsearch",action="store_true",dest="dmsearch",                         help="search DMs with prepfold")
    pars.add_argument("--psearch",action="store_true",dest="psearch",                           help="search period with prepfold")
    pars.add_argument("--nonew",action="store_true",dest="nonew",                               help="do not search for .fits files")
    pars.add_argument("--ignorechans",action="store",dest="igchn",                              help="file with channels to ignore")
    pars.add_argument("--pos",action="store",                                                   help="ra and dec ('hms dms' or 'hms:dms')")
    args = vars(pars.parse_args())

    if args['pos'] != None:
        if len(args["pos"].split()) > 1:
            pos=SkyCoord(args['pos'],unit=('hourangle', 'deg'))
        elif len(args['pos'].split(':')) > 1:
            pos=SkyCoord(args['pos'].split(':')[0],args['pos'].split(':')[1],unit=('hourangle','deg'))
    args["ra_psr"] = pos.ra.value
    args["dec_psr"] = pos.dec.value
    return args

args = proc_args()

p0_psr          = args["period"] 	# (N) period of command line-input pulsar (will be ignored if file is provided)
pd_psr          = args["pdot"]		# (N) period derivative of command line-input pulsar (will be ignored if file is provided)
dm_psr          = args["dm"]		# (N) dm of command line-input pulsar (will be ignored if file is provided)
ra_psr		= args["ra_psr"]	# (N) ra of command line-input pulsar (will be ignored if file is provided)
dec_psr		= args["dec_psr"]	# (N) dec of command line-input pulsar (will be ignored if file is provided)
beam_psr        = args["beam"]		# (N/F) specific beam to check for command line-input pulsar (will be ignored if file is provided) / file with columns JNAME BEAM[comma-sep] from which to read beams 
file_psr        = args["file"]		# (F) input file from which pulsars will be read (name rajd decjd p0 p1 dm); preferably in astropy table format
nbin            = args["nbin"]		# (N) number of bins to be used in prepfold pulse profile (otherwise, will choose based on period)
nsub            = args["nsub"]		# (N) number of subbands to used in prepfold frequency domain (defaults to 128)
npart           = args["npart"]		# (N) number of subintegrations in prepfold (defaults to 64)
ang_max         = args["angle"]		# (N) maximum angular offset between beams and pulsar (defaults to 0.5 degrees)
snr_min         = args["snr"]		# (N) threshold for detection
closest         = args["closest"]	# (B) only use beam closest to pulsar
data_dir        = args["dat"]		# (D) directory where .fits files live; defaults to /lustre/cv/projects/GBNCC/
point           = args["point"]		# (F) input file from which survey beam positions are stored; preferably in astropy table format; beam rajd decjd
rfi_fil         = args["rfi"]		# (B) use additional rfi-mitigation flags (zapping 2470:3270 channels after MJD 56710, -noscales, -noweights, adjusting some rfifind settings
par_file        = args["par"]		# (F/D) par file or directory where par file is stored
rewrite         = args["rewrite"]	# (B) boolean for deleting previous results and redoing them
outrt           = args["outputroot"]    # (S) root name for summary file
work_dir        = args["out"]		# (D) directory where _temp files will be created and processing will occur
dmsearch        = args["dmsearch"]	# (B) turn on DM searching
psearch         = args["psearch"]	# (B) turn on period searching
nonew           = args["nonew"]		# (B) don't look for .fits files (use after processing at least once without this flag)
igchn           = args["igchn"]		# (F) file with channels to ignore for pulsars (list does not have to match input file; JNAME CHANS[comma-sep])
use_comp        = 'zuul'		# (S) in the future, can be made an option to look in different places


# Reference tables
if point.split('.')[-1]=='fits':
    point=Table.read(point)
else:
    point=Table.read(point,format='ascii')
if len(point['ra'][0].split(':'))>1:
    surv=SkyCoord(point['ra'],point['dec'],unit=('hourangle','deg'))
else:
    surv=SkyCoord(point['ra'],point['dec'],unit='deg')
beam_rejects=[]
beam_prefix='GBNCC'
#beam_prefix='GBT820-'
cent_freq=350.
#cent_freq=820.

if str(file_psr)!='None'
    if file_psr.split('.')[-1]=='fits':
	psr_list=Table.read(file_psr)
    else:
	psr_list=Table.read(file_psr,format='ascii')
else:
    ps=SkyCoord(ra_psr,dec_psr,unit='deg').split()
    name='J'+str(ps[0].split('m')[0][0:2])+str(ps[0].split('m')[0][3:])+str(ps[1].split('m')[0][0:3])+str(ps[1].split('m')[0][4:])+'_cand'
    try:
	psr_list=Table(names=('name','rajd','decjd','p0','p1','dm'),dtype=('S10','float','float','float','str','float'))
	psr_list.add_row([name,ra_psr,dec_psr,p0_psr,p1_psr,dm_psr])
    except:
	print "problem with input; check data types"
	exit()

cond=dot(dot(psr_list['decjd']>=surv.dec.value.min()-ang_max,psr_list['decjd']<=surv.dec.value.max()+ang_max),dot(psr_list['rajd']>=surv.ra.value.min()-ang_max,psr_list['rajd']<=surv.ra.value.max()+ang_max))
print "%i/%i valid pulsars in survey area" %(sum(cond),len(cond)
psr_list=psr_list[cond]

if "outputroot" in args:
    summary_file = open("%s%s_%s_%dpulsars.txt" %(work_dir,outrt,use_comp,sum(cond)), 'w')
else:
    summary_file = open("%sbeamcheck_%s_%dpulsars.txt" %(work_dir,use_comp,sum(cond)), 'w')
summary_file.write("#@ Using S/N > %s for detection threshold on %s\n" %(snr_min,use_comp))
if str(beam_psr) == "None" and file_psr is not None:
    summary_file.write("#@ Using beams within %.1f deg of pulsars in %s \n" %(ang_max,file_psr))
elif str(beam_psr) != "None" and file_psr is not None:
    summary_file.write("#@ Using beam #%s to search for pulsars from %s \n" %(beam_psr,file_psr))
elif str(beam_psr) != "None" and file_psr is None:
    summary_file.write("#@ Using beam #%s to search for pulsar %s \n" %(beam_psr,pulsar.name))
elif str(beam_psr) == "None" and file_psr is None:
    summary_file.write("#@ Using beams within %.1f deg of pulsar %s \n" %(ang_max,pulsar.name))
summary_file.write("#@ Flux adjusted for offset from beam center\n")
summary_file.write("#@ Flux estimated from measured S/N, uncertainty estimated from dependencies\n")
summary_file.write("#@ File generated %s (ET)\n" %strftime("%Y-%m-%d %H:%M:%S"))
if dmsearch:
    summary_file.write('#@ DM Searching - DMs is resulting DM from search\n')
    summary_file.write("#@  PSRJ RA Dec P0 DM DMs Dist MJD BW Temp S/N Beam Status\n")
    summary_file.write("#@ -- (deg) (deg) (s) (p/cc) (p/cc) (deg) (day) (MHz) (K) -- -- --\n")
else:
    summary_file.write("#@ PSRJ RA Dec P0 DM Dist MJD BW Temp S/N Beam Status\n")
    summary_file.write("#@ -- (deg) (deg) (s) (p/cc) (deg) (day) (MHz) (K) -- -- --\n")

# Tracking number for final output
name_detect = []
num_detect = 0
num_unobs = 0

mult=0
for l in psr_list:
    psr=Pulsar([l['name'],l['rajd'],l['decjd'],l['p0'],l['p1'],l['dm']])
    if sum(psr_table['name']==psr.name)>1:
	print "multiple pulsars with name %s - appending _[index] to names"
	psr.name=psr.name+'_%i' %mult
	mult+=1
    pos=SkyCoord(psr.ra,psr.dec,unit='deg')
    if str(beam_psr)=='None':
	if not closest:
	    for b in point[pos.separation(surv).value<=ang_max]['beam']:
		if b not in beam_rejects:
		    psr.beams.append(b)
	else:
	    psr.beams.append(Beam(point[pos.separation(surv)==pos.separation(surv).min()]['beam']))
    else:
	if os.path.isfile(beam_psr):
	    for b in sproc.Popen("grep %s %s" %(psr.name.split('_')[0],beam_psr),shell=True,stdout=sproc.PIPE).communicate()[0].split()[1:]
		psr.beams.append(Beam(b))
	else:
	    psr.beams.append(Beam(beam_psr))
    print psr.name
    psr.beam_numbers()
    if str(par_file)!='None':
        if os.path.isfile(par_file):
	    psr.par=par_file
        elif os.path.isdir(par_file):
	    psr.par=glob('%s/%s*par' %(par_file,psr.name))[0]
        else:
	    print "unable to find %s" %par_file
    if psr.par='':
        if len(glob('/users/rspiewak/pulsars/%s*par' %psr.name))>0:
	    psr.par=glob('/users/rspiewak/pulsars/%s*par' %psr.name)[0]
	if len(glob('/lustre/cv/projects/GBNCC/amcewen/pars/%s*par' %psr.name))>0:
	    psr.par=glob('/lustre/cv/projects/GBNCC/amcewen/pars/%s*par' %psr.name)[0]
    if psr.par='':
	print "cannot find par file. using input parameters for folding"
    os.chdir(work_dir)
    if not os.path.isdir('%s_temp' %psr.name):
	os.mkdir('%s_temp' %psr.name)
    elif rewrite:
	sproc.call('rm %s_temp/*' %psr.name,shell=True)
    os.chdir('%s_temp' %psr.name) 
    if psr.par!='':
	sproc.call('ln -sf %s .' %psr.par,shell=True)
    for beam_cand in psr.beams:
	if os.getcwd() != "%s%s_temp" %(work_dir,psr.name):
	    os.chdir("%s%s_temp" %(work_dir,psr.name))
	beam_cand.ang_off(psr.ra,psr.dec)
	if len(glob('*%s_*fits %beam_cand.name))==0:
	    fits_list=[]
	    if not nonew: 
		fits_list=glob('%s2*/gu*%s_*fits' %(data_dir,beam_cand.name))
	    if nonew or len(fits_list)==0:
		print "Cannot find beam %s for psr %s" %(beam_cand.num,psr.name)
		if dmsearch:
		    summary_file.write("%s %s %s %s %s -- %s -- -- %s %s %s NOTFOUND" %(psr.name,psr.ra,psr.dec,psr.p0,psr.dm,beam_cand.off,psr.temp,beam_cand.num))
		else:
		    summary_file.write("%s %s %s %s %s %s -- -- %s %s %s NOTFOUND" %(psr.name,psr.ra,psr.dec,psr.p0,psr.dm,beam_cand.off,psr.temp,beam_cand.num))
		num_unobs+=1
		os.chdir(work_dir)
		continue
	    for i in fits_list:
		sproc.call('ln -sf %s .' %i,shell=True)
	for f in glob('*%s_*fits %beam_cand.name):
	    beam_cand.add_fits(f.split('/')[-1])
	rfi_std = "-time 1 -timesig 2 -freqsig 3"
	for fits,mjd in zip(beam_cand.fits,beam_cand.mjd):
	    if os.getcwd().split('/')[-1]!="%s_temp" %psr.name:
		os.chdir("%s_temp" %psr.name)
	    if '2bit' in fits and rfi_fil:
		raw_opt = "-noscales -noweights"
		rfi_std = "-time 1 -freqsig 3"
	    else: 
		raw_opt = ""
	    if mjd>56710 and rfi_fil:
		rfi_opt = "-zapchan 2470:3270"
	    else:
		rfi_opt = ""
	    if len(glob('*%s*%s*rfifind.mask' %(mjd,beam_cand.num)))==0 or len(glob(*%s*%s*rfifind.stats' %(mjd,beam_cand.num)))==0:
		if len(glob('%smask_files/*%s*%s*rfifind.mask' %(data_dir,mjd,beam_cand.num)))==0 or len(glob('%smask_files/*%s*%s*rfifind.stats' %(data_dir,mjd,beam_cand.num)))==0:
		     

