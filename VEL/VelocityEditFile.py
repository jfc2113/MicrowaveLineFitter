import numpy as np, pylab as pl, csv

#####################################################################################################
##############################  This is the one you edit!              ##############################
#####################################################################################################


#####################################################################################################
#####################  This script will conver the output of FullCode.py to a concise csv file with 
#####################  line IDs and gaussian fits in velocity space.
#####################  Before running, re-fit any lines in ALLFITS with bin = 5, and fix any other 
#####################  fits as desired.  Also, prepare 1-2 csv files containing prioritized line IDs.  
#####################  As it is currently written, the code inspects assigned recombination lines (priority 1)
#####################  for kinematic consistency first.  The script does not require a csv file input 
#####################  for the recombination lines at this stage however.  Then (priority 2) it looks for 
#####################  lines that you specify as expecting to be strong, in a file "strongLines.csv".
#####################  Third, it handles  foreground absorption by very strong lines held in 'SAClines.csv'
#####################  (SAC stands for spiral arm cloud).  Fourth, it assigns lines detected in ALLFITS.csv.  
#####################  Finally, it handles unidentified transitions.
#####################  Users can specify either 1 or 2 csv files to establish the priorities as described below.

def REGIONS(go):
	return ['K6','L','K4']		##EDIT THIS LINE TO INDICATE THE REGIONS YOU WILL INVESTIGATE.  FOR EACH REGION (r), YOU WILL NEED A FILE CALLED "../ALLFITS_(r).csv"

## Specify the input csv files including the output of FullCode.py and any files containing line frequencies.  


def ALLFITSin(r):
	return '../ALLFITS_%s.csv' % r, 'GHz'			## REQUIRED  # Edit the name of the file as needed.  Include the units that the line rest frequencies are listed in.  (These are in Column J of the standard ALLFITS file.)

def velshift(r):	
	return 64	#vLSR in km/s		## ARE YOUR DATA ALREADY VELOCITY SHIFTED? IF NOT, CHANGE TO 0.  If towards different regions, different vels, make this a function of r.

def strongAndSAC(r):
	## Optional:  Point to a csv file containing strong molecular lines expected in the source, lines that may appear in foreground absorption (in translucent & diffuse clouds), and the unit in which the frequencies are listed.  If you do not want to prepare these files, change the next line to "return []"
	strong = 'strongLines.csv'
	foreground='SAClines.csv'
	units = 'GHz'
	return [strong, foreground, units]
	## Set any of these to 'dum' or any other placeholder if you dont want to use these csvfiles.
	## Another example:  return ['mySpecialLines.csv','TranslucentChem.csv','MHz']
	## This could also be made to vary by the source:  if r == 'OrionKL': return ['hotCores.csv','dum','GHz']; if r == 'OrionBar': return ['PDRs.csv','dum','GHz']; if r == 'W49N': return ['hotCores.csv','TranslucentChem.csv','GHz']

def goodBins(r):
	return [0,1,2,3,4]	## These are the bins that have valid fits and will be output into the velocity_ALL file.  If you only want to include very good lines, include [0,1,2].  bin=4 is what I use to indicate I've manually re-fit a line

##########  PROVIDE STRONG MASER FREQUENCIES IF YOU HAVE RINGING OR STRANGE LINE SHAPES AROUND THESE LINES.  ###########


def RelevantMaserFreqs(r): 
	masefreqs = np.asarray([36169.29])		#in MHz, the super strong methanol maser frequency.  Because of ringing & due to the narrow line widths, these lines have many wings extending further away than weaker, normal transitions.  Strong masing lines should be treated separately in selecting 'wing' components.
	velocityWingRange=50			# this is the velocity range around a masing transition in which we will readily assign wing components
	return masefreqs, velocityWingRange	


##########  EDIT THESE BASED ON THE SOURCES YOU ARE OBSERVING  ###########

def SACvRange(r):		
## This outputs a velocity range that is reasonable for foreground absorption by translucent or diffuse cloud material.  Only applies to the strongest lines, which you should specify with the parameter "SACfile" above.
	if r =='L':			Vrange=[-120,50]
	if r =='LMH':			Vrange=[-120,55]
	if r == 'K6' or r == 'N':	Vrange=[-120,58]
	if r =='K4': 			Vrange=[-20,60]
	if r =='M':			Vrange=[-120,40]
	return Vrange

def strongVrange(r):		
## This outputs a velocity range that is reasonable for strong lines detected in the source.  Only applies to lines specified with the parameter "strongLines" above.
	if r =='L':		
		Vrange=[50,80]
		vmain = 55
	if r =='LMH':
		Vrange=[55,85]
		vmain = 64
	if r == 'K6' or r == 'N': 
		Vrange =[58,85]
		vmain = 64
	if r =='K4':
		Vrange =[58,88]
		vmain = 62
	if r == 'h':
		Vrange =[60,88]
		vmain = 73
	if r == 'M':
		Vrange =[40,70]
		vmain = 50
	return Vrange, vmain

def bestcheck(vel,region,linetype='dum',height=0,velcr='dum',velwid='dum'):
## This specifies the "best" velocity range.  You could use similar parameters to those from strongVrange. 
	if velcr!='dum' and velwid!='dum':
		if abs(vel - velcr) < abs(velwid):
			return 'True'

	if linetype =='Recomb' or linetype == 'recomb':
		if region =='K6':
			if abs(vel-75 < 18):
				return 'True'

		else:
			if abs(vel-70) < 15:
				return 'True'
			else:
				return 'False'
		
	if region =='L':
		if linetype =='mol':		
			return 'True'
		if linetype == 'strong': 
			if abs(vel - 55.) <5:
				return 'True'
			if abs(vel - 75.) <5:
				return 'True'
			else:
				return 'False'
		if linetype =='maser':
			if height > 500.:
				return 'True' 
		if linetype == 'SAC':  
			if abs(vel + 50) < 70:  
				return 'True'
			else:
				return 'False'
				
	if region =='LMH':
		if linetype =='mol':	
			if abs(vel - 60.) <4:
				return 'True'
			else:
				return 'False'
		if linetype == 'strong': 
			if abs(vel - 64.) <3 or  abs(vel - 82.) <3:
				return 'True'
			else:
				return 'False'
		if linetype =='maser':
			if height > 500.:
				return 'True' 
			else:
				return 'False'
		if linetype == 'SAC':  
			if abs(vel + 50) < 70:  
				return 'True'
			else:
				return 'False'
	if region =='h':
		if linetype =='mol':	
			if abs(vel - 71.) <4:
				return 'True'
			else:
				return 'False'
		if linetype == 'strong': 
			if abs(vel - 64.) <3 or  abs(vel - 82.) <3:
				return 'True'
			else:
				return 'False'
		if linetype =='maser':
			if height > 500.:
				return 'True' 
			else:
				return 'False'
		if linetype == 'SAC':  
			if abs(vel + 50) < 70:  
				return 'True'
			else:
				return 'False'
	if region =='N' or region == 'K6':
		if linetype =='mol':	
			if abs(vel - 64.) <4:
				return 'True'
			else:
				return 'False'
		if linetype == 'strong': 
			if abs(vel - 64.) <4 or  abs(vel - 82.) <3:
				return 'True'
			else:
				return 'False'
		if linetype =='maser':
			if height > 500.:
				return 'True' 
			else:
				return 'False'
		if linetype == 'SAC':  
			if abs(vel + 50) < 70:  
				return 'True'
			else:
				return 'False'
	if region =='M':
		if linetype =='mol':	
			if abs(vel - 70.) <10:
				return 'True'
			else:
				return 'False'
		if linetype == 'strong': 
			if abs(vel - 62.) <6 or  abs(vel - 78.) <6:
				return 'True'
			else:
				return 'False'
		if linetype =='maser':
			if height > 500.:
				return 'True' 
			else:
				return 'False'
		if linetype == 'SAC':  
			if abs(vel + 50) < 70:  
				return 'True'
			else:
				return 'False'
	if region =='K4':
		if linetype =='mol':	
			if abs(vel - 63.) <4:
				return 'True'
			else:
				return 'False'
		if linetype == 'SAC':  
			if abs(vel - 10) < 20:  
				return 'True'
			else:
				return 'False'
		if linetype == 'strong': 
			if abs(vel - 60.) <6.5 or  abs(vel - 82.) <2:
				return 'True'
			else:
				return 'False'
		if linetype =='maser':
			if height > 500.:
				return 'True' 
			else:
				return 'False'
	else:
		if linetype != 'wi' and linetype != 'wing' and linetype != 'win':
			print "gotta give bestcheck more inputs to go on."
			print "vel,region,linetype,height,velcr,velwid =", vel,region,linetype,height,velcr,velwid
		return 'False'

