#################  FullCode.py  ################
##   This is the automated line fitter & identifier. 
##   All csv files used have the following settings:  Character set = Western Europe (Windows-1252/WinLatin 1); Language = Default-English (USA); separated by = :



# IMPORT DEPENDENCIES
import numpy as np, pylab as pl
from math import *
import LineFind as LF, Gfit as Gfit, LineID as lID, FTinterp as fti
c=2.9979E5		## Speed of light (km/s)
## LineFind is used to grab segments of the spectrum that contain lines.  Gfit is used to do the line fitting.  LineID compares to specified csv files that contain the output of a line catalog.  FTinterp performs interpolation by Fourier-domain zero padding. 

##############################################################################################
###############################   EDIT THIS SECTION AS NEEDED   ##############################
##############################################################################################
thresh = 3.5  		## This many sigma is the threshold for line detection.
wayout=False  		## If wayout == 'True', after determining the regions with lines, the code will plot the full spectrum and the line regions without fitting anything.  Can use as a check before doing line fits, or if you see something odd in the fits.
writeOut = 'True'  	## Do you want to write out a csv file containing the line fit parameters?
wantLineIDs=True	## Do you want to do the line identification?  Doing so requires that you have an input csv file (compFile1, below) that contains line information.  
plotLineIDs=True	## Do you want to plot the frequencies of lines from catalogs that correspond to your detected transition? This also requires that you have an input csv file (compFile1, below) that contains line information. Won't happen if wantLineIDs == False. 
GoutPlotQ = False	## If GoutPlotQ == True, separate plots will be made for each line that is fit.  For doing the full spectrum, GoutPlotQ = True will generate MANY plots!!Generally will want GoutPlotQ = False, which will generate a single plot of the full spectrum with the full fit overlaid.
GQuiet = True		## If GQuiet == False, will print out multiple updates while running the Gaussian Fitter (Gfit.gauFit).
ExtraQuiet = True	## If ExtraQuiet == False, will print out a few updates during the loop below.
plotLID_Height='dum'	## If you plot line IDs, this specifies the height you want to plot.  
LenGauArs = 100000	## This is the length of the array that holds the Gaussian fit to the spectrum.  The "resolution" of the line-fit spectrum is (Freq_Max-Freq_Min)/LenGauArs.  Increase for very broadband work, such that LenGauArs > 6 * Total_Bandwidth/Line_Width  
wantInterp = True	## If this is true, the Gaussian fitter will fit interpolated data instead of the raw data.  Use if your lines are poorly sampled.  Set to False if your data is <2x Nyquist sampled.
startCRF,endCRF,sepVal = 30750,49150,1850	## If you have a spectrum consisting of multiple sections (multiple spectral windows or images), startCRF is the the center frequency of the first section, endCRF is the center frequency of the last section, and sepVal is the separation between section (all in MHz).
tooBroadvAbs = 35	## Set a width that is too broad to be reasonable for an absorption line.  In km/s.  Default is set to 35km/s.
tooBroadV='dum'		## Set a width that is too broad to be reasonable for an emission line.   In km/s.  Default is set to 70km/s.
broadV='dum'		## Set a width that is broader than expected, but not utterly unreasonable for an emission line.  In km/s. Default is set to 40 km/s
tooNarrowV = 'dum'	## Set a width that is too narrow to be reasonable for an emission line.  In km/s.  Default is set to 2.2km/s.  If your lines are narrow, need to change!
compSpread=20		## This is the typical velocity separation between velocity components in km/s.  Default is set to 20km/s.
compWid=15		## This is the typical width of velocity components.  Default is set to 15km/s.
dynRange=100		## This is the dynamic range of your measurements.  For example, in the immediate vicinity of a very strong masing transition, the dynamic range of the data set is ~100, even though the ratio of the masing transition strength to the weakest detectable line in the image (but at a slightly separated frequency), is much larger.
velShift=64		## If your data is already shifted by some velocity, put this here.  In this case, the spectrum is already shifted to +64 km/s.

## THE SPECTRUM WILL BE COMPARED WITH LINE PARAMETERS FROM A CATALOGUE.  SET THE LINE PARAMETER FILES HERE.  A SPLATALOGUE OUTPUT CSV FILE (with a ':' delimeter) REQUIRES MINOR FORMATTING TO BE INPUT TO THIS CODE.
## Note:  Keep recombination lines and molecular lines in seperate files.
regions=['K6', 'L', 'K4']			#Set the names of the regions to be fit.
compFile1 = 'recombLines.csv'		#Contains recombination lines used for line-IDs.  If you expect no recomb lines, Use compFile1 = 'molecularLines.csv', and have this file contain possible line identifications.  See 'recombLines.csv' for formatting requirements.  Do not make one file with both recomb and molecular lines -- keep these seperate.
compFile2 = 'molecularLines.csv'	#Contains molecular lines used for line-IDs.  The code is set to do these separately from the recomb lines because the recomb lines are very reliably detected = 'molecularLines.csv', and have this file contain possible line identifications.  Set compFile2 = 'None' or compFile2='dum' if you are only comparing to one file.  See 'recomb.csv' or 'molecularLines.csv' for their respective formatting requirements.  
compFileUnits = 'GHz'			#The frequencies in recomb.csv & molecularLines.csv are in units of GHz.  Change if yours are in units of MHz.  Have all compFiles in the same units.

### The data input is defined here.
### The spectrum should be an ascii file with a frequency column and an intensity column.  The frequency column should be in MHz.  
def InSpec(r):
	return LF.spect('%s_fullSpec.txt' %(r),units='MHz')  	#Grabs a spectrum named 'K6_ALL.txt' if K6 is the current region, r. Specify if your spectrum is in MHz or GHz.
yunits = 'Specific Flux Density (mJy arcsec$^{-2}$)'			##Specify the units of the intensity units in the file in InSpec.  

### The csv file output is defined here.  
def OF(r):
	return "ALLFITS_%s.csv" % r				#Writes a csv file named 'ALLFITS_K6.txt' if K6 is the current region, r. This file will be built to contain all of the line fit parameters and line identifications.

def LineFind_velParams(r):
	if r == 'K6':
		recombVel=75	# Center velocity for selecting recombination lines.
		recombdelVel=25	# +- velocity range for identifying recombination lines.  All recomb lines within [25] km/s will be noted in the output lineID file.
		molVel=73	# Center velocity for selecting molecular lines
		moldelVel=17	# +- velocity range for identifying molecular lines.   
		molBest=64	# the best velocity for molecular lines.
	if r =='L':
		recombVel=75
		recombdelVel=15
		molVel=64
		moldelVel=17
		molBest=55
	if r =='K4':
		recombVel=80
		recombdelVel=25
		molVel=72
		moldelVel=17
		molBest=62
	return recombVel,recombdelVel,molVel,moldelVel,molBest

##############################################################################################
#########################   SHOULD NOT NEED EDITS BELOW THIS POINT.  #########################
##############################################################################################

## This function returns interpolated data if wantInterp is set to True.
def interps(sFreq,sFlux,wantI):
	if wantI == True or wantI == 'True':
		iFreq,iFlux=fti.padding(sFreq,sFlux,plotQ='False')  ## Performing interpolation by Fourier-domain zero padding.
	if wantI == False or wantI == 'False':
		iFreq,iFlux=sFreq,sFlux
	return iFreq,iFlux

## This function fills a set of arrays that are output to the output csv file.
def listfill(lineAr,hAr,crAr,widAr,FMAr,croffAr,binAr,lineC,Goutp,freqAr,GT,bigFreq):
	if len(Goutp)==6:
		lineAr.append(lineCounts)
		hAr.append(Goutp[1])
		FMAr.append(Goutp[4])
		croffAr.append(Goutp[5])
		widAr.append(Goutp[3])
		crAr.append(Goutp[2])
		binAr.append(Goutp[0])
		GG=Gfit.gauss_fit1(Goutp[1:4],freqAr)
		GT+=Gfit.gauss_fit1(Goutp[1:4],bigFreq)
	if len(Goutp)==8:
		lineAr.extend([lineC,lineC])
		hAr.extend([Goutp[1],Goutp[4]])
		crAr.extend([Goutp[2],Goutp[5]])
		widAr.extend([Goutp[3],Goutp[6]])
		FMAr.extend([Goutp[7],Goutp[7]])
		binAr.extend([Goutp[0],Goutp[0]])
		croffAr.extend([1000,1000])
		GG=Gfit.gauss_fit2(Goutp[1:7],freqAr)
		GT+=Gfit.gauss_fit2(Goutp[1:7],bigFreq)
	return lineAr,hAr,crAr,widAr,FMAr,croffAr,binAr,GG,GT

#####################################################################################################################
####################   Now run through a loop for each spectrum.   ##################################################
#####################################################################################################################
for r in regions:
	outfile = OF(r)
	MoreGs=0
	fig = pl.figure(figsize=(12,6))
	ax = fig.add_subplot(111)
	lineCounts=0
	f,fd = InSpec(r)   					## The input spectrum.
	fRes = np.median(f[1:]-f[:len(f)-1])			## The channel resolution.	
	if startCRF<f.min():					## If so, the startCRF, endCRF, and sep values were probably not set.  Defaults to treating the entire spectrum as having uniform noise.
		startCRF = f.mean()
		endCRF = f.mean()+(f.max()-f.mean())/2
		sepVal = f.max()-f.mean()
		print "****  startCRF, endCRF, and sepVal were not set properly.  Treating the entire spectrum as if it has uniform noise  ****"
	freq,flux,regRMS=LF.LineFind(f,fd,plotQ='False',rmsthresh1=thresh,startCR=startCRF,endCR=endCRF,sep=sepVal)	## Determines where lines are.  The freq & flux arrays ONLY contain segments of the spectrum with lines.
	if wayout ==True or wayout == 'True':	## If wayout == 'True', after determining the regions with lines, the code will plot the full spectrum and the line regions without fitting anything.  Can use as a check before doing line fits or if you see something odd in the fits.
	   print "wayout == True.  Will show the regions determined to detect lines and exit instead of fitting the spectrum."
	   pl.plot(f,fd,'c',linewidth=5,label='Original Spectrum')
	   starting = 0
	   freqDif=freq[1:]-freq[:len(freq)-1]
	   dumdum = 0
	   for iiii in range((len(freqDif))):
		if abs(freqDif[iiii])>1.9:
			ending = iiii+1
			if ending-starting > 3:
				smallFreq=freq[starting:ending]	
				smallFlux=flux[starting:ending]
				if dumdum ==0:
					pl.plot(smallFreq,smallFlux,color='y',linewidth=5, label='Segments Containing Lines')
					if wantInterp == True:
						interpFreq,interpFlux=fti.padding(smallFreq,smallFlux,plotQ='False')				
						pl.plot(interpFreq,interpFlux,color='b',linewidth=1, label= 'Interpolated Segments')
					dumdum=1
				else:
					pl.plot(smallFreq,smallFlux,color='y',linewidth=5)
					if wantInterp == True:
						interpFreq,interpFlux=fti.padding(smallFreq,smallFlux,plotQ='False')				
						pl.plot(interpFreq,interpFlux,color='b',linewidth=1)
			starting = ending

  	   pl.xlabel('Frequency (MHz)')
	   pl.ylabel(yunits)
	   pl.legend(prop={'size':8})
	   pl.show()
	   break

	GauTotal = np.zeros((LenGauArs))
	GauXaxis = np.linspace(f.min(),f.max(),LenGauArs)
	regRMS0=regRMS
	starting = 0
	freqDif=freq[1:]-freq[:len(freq)-1]					
	hs, FMs, croffs, crs, wids, binCodes,lineIDs=[],[],[],[],[],[],[]	## These will hold the Gaussian line fit parameters.
	for iiii in range((len(freqDif))):		## Run through the array containing segments of the region with lines.
		if abs(freqDif[iiii])>1.9*fRes:		## If the frequency axis indicates a break significantly larger than a channel, it enters this loop to fit the data in one segment of the spectrum.
			ending = iiii+1
			if ending-starting > 3:		## It won't try to fit segments that are smaller than 3 channels wide.  Should only matter at the edges of the spectrum.	
				if GQuiet == False or GQuiet == 'False':
					print "iiii= ", iiii
				smallFreq, smallFlux = freq[starting:ending], flux[starting:ending]	## Grabbing a single segment of the spectrum.
				interpFreq,interpFlux=interps(smallFreq,smallFlux,wantInterp)  ## Interpolating segment of spectrum if wantInterp == True.
				if ExtraQuiet == 'False' or ExtraQuiet == False: 
					print " \n \n \n"
					print "In FullCode:  Passing Freq_min = %.2f to Freq_max =%.2f to Gfit" % (interpFreq.min(),interpFreq.max())
				Gout=Gfit.GauFit(interpFreq,interpFlux,regRMS,1,[0,0,0],quiet=GQuiet, plotQ=GoutPlotQ,tooBroadvABS=tooBroadvAbs,Broadv=broadV,tooBroadv=tooBroadV,tooNarrowv=tooNarrowV,compSpread=compSpread,compWid=compWid,DynRange=dynRange,firstImFreq=startCRF,imFreqRange=sepVal)	## Broad absorption is ~5MHz wide.
				if ExtraQuiet == False or ExtraQuiet == 'False':
					print "Gout[0]",Gout[0]
				lineIDs,hs,crs,wids,FMs,croffs,binCodes,GGG,GauTotal = listfill(lineIDs,hs,crs,wids,FMs,croffs,binCodes,lineCounts,Gout,smallFreq,GauTotal,GauXaxis)
				lineCounts+=1			
				for dumi in np.arange(2,20,1):
					smallFlux=smallFlux-GGG
					ival = int((Gout[2]-29825.)/1850)
					if ival+1 > len(regRMS) or ival < 0:
					    if ival <0:
						regSpecRMS = regRMS[0]
					    if ival+1 > len(regRMS):
						regSpecRMS = regRMS[ival-1]
					else:
						regSpecRMS = regRMS[ival]
					if abs(Gfit.maxabs(smallFlux))< thresh*regSpecRMS:
						GGG = np.zeros((len(GGG)))
						break	
					if ExtraQuiet == False or ExtraQuiet == 'False':
						print "Fitting G in iteration %i" % dumi
					interpFreq,interpFlux=fti.padding(smallFreq,smallFlux,plotQ='False')
					Gout=Gfit.GauFit(interpFreq,interpFlux,regRMS,1,[0,0,0],quiet=GQuiet, plotQ=GoutPlotQ,tooBroadvABS=tooBroadvAbs,Broadv=broadV,tooBroadv=tooBroadV,tooNarrowv=tooNarrowV,compSpread=compSpread,compWid=compWid,DynRange=dynRange,firstImFreq=startCRF,imFreqRange=sepVal)	
					if ExtraQuiet == False or ExtraQuiet == 'False':
						print "Gout[0]",Gout[0]
					lineIDs,hs,crs,wids,FMs,croffs,binCodes,GGG,GauTotal = listfill(lineIDs,hs,crs,wids,FMs,croffs,binCodes,lineCounts,Gout,smallFreq,GauTotal,GauXaxis)
					lineCounts+=1
			starting = ending

	smallFreq=freq[starting:len(freq)]	
	smallFlux=flux[starting:len(freq)]
	interpFreq,interpFlux=fti.padding(smallFreq,smallFlux,plotQ='False')
	if ExtraQuiet == 'False' or ExtraQuiet == False: 
		print "  "
		print "In FullCode:  Passing interpFreq.min() = %.2f to interpFreq.max() =%.2f to Gfit" % (interpFreq.min(),interpFreq.max())
	Gout=Gfit.GauFit(interpFreq,interpFlux,regRMS,1,[0,0,0],quiet=GQuiet, plotQ=GoutPlotQ,tooBroadvABS=tooBroadvAbs,Broadv=broadV,tooBroadv=tooBroadV,tooNarrowv=tooNarrowV,compSpread=compSpread,compWid=compWid,DynRange=dynRange,firstImFreq=startCRF,imFreqRange=sepVal)	
	lineIDs,hs,crs,wids,FMs,croffs,binCodes,GGG,GauTotal = listfill(lineIDs,hs,crs,wids,FMs,croffs,binCodes,lineCounts,Gout,smallFreq,GauTotal,GauXaxis)

	for dumi in np.arange(2,20,1):
		smallFlux=smallFlux-GGG
		ival = int((Gout[2]-29825.)/1850)
		if ival+1 > len(regRMS) or ival < 0:
		    if ival <0:
			regSpecRMS = regRMS[0]
		    if ival+1 > len(regRMS):
			regSpecRMS = regRMS[ival-1]
		else:
			regSpecRMS = regRMS[ival]
		if abs(Gfit.maxabs(smallFlux))< 3*regSpecRMS:
			GGG = np.zeros((len(GGG)))
			break	
		if ExtraQuiet == False or ExtraQuiet == 'False':
			print "Fitting G in iteration %i" % dumi
		interpFreq,interpFlux=fti.padding(smallFreq,smallFlux,plotQ='False')
		Gout=Gfit.GauFit(interpFreq,interpFlux,regRMS,1,[0,0,0],quiet=GQuiet, plotQ=GoutPlotQ,tooBroadvABS=tooBroadvAbs,Broadv=broadV,tooBroadv=tooBroadV,tooNarrowv=tooNarrowV,compSpread=compSpread,compWid=compWid,DynRange=dynRange,firstImFreq=startCRF,imFreqRange=sepVal)	
		if ExtraQuiet == False or ExtraQuiet == 'False':
			print "Gout",Gout
		lineIDs,hs,crs,wids,FMs,croffs,binCodes,GGG,GauTotal = listfill(lineIDs,hs,crs,wids,FMs,croffs,binCodes,lineCounts,Gout,smallFreq,GauTotal,GauXaxis)
		lineCounts+=1
		
	regionMaster=zip(lineIDs,crs,croffs,hs,wids,FMs,binCodes)
	crs=np.asarray(crs)
	if ExtraQuiet == False or ExtraQuiet == 'False':
		print "\n################################\nDone with Gaussian Fitting."

####################   Now output and plot line identifications if plotLineIDs is True.   ###############################################
#####################################################################################################################
	if writeOut == True or writeOut == 'True' or plotLineIDs == 'True' or plotLineIDs==True or wantLineIDs==True or wantLineIDs=='True':	

		if ExtraQuiet == False or ExtraQuiet == 'False':
			print "Now beginning the Line-IDs."
		rV,rdV,molV,moldV,molBest = LineFind_velParams(r)
		if wantLineIDs == True or wantLineIDs=='True':
			lID.lineCompareAll(outfile,regionMaster,infile1 = compFile1,infile2 = compFile2,plotQ=plotLineIDs,writeout=writeOut,velshift=velShift,velRecomb=rV, dV_recomb=rdV, velMol = molV,dV_mol=moldV,best=molBest,plotlineHeight=plotLID_Height,inFileUnits=compFileUnits)
		else:
			lID.lineCompareAll(outfile,regionMaster,plotQ=plotLineIDs,writeout=writeOut,velshift=velShift,velRecomb=rV, dV_recomb=rdV, velMol = molV,dV_mol=moldV,best=molBest,inFileUnits=compFileUnits,plotlineHeight=plotLID_Height)
			
		print "should have written outfile = %s out" % outfile


	for icrfreqi in range(len(regRMS)):
		crf = startCRF+sepVal*icrfreqi
		rmsNow= regRMS[icrfreqi]
		print "for image at crFreq = ",crf, "RMS = ",rmsNow 
	pl.plot(f,fd,color = 'c')#,linewidth=4)
	interpFreq,interpFlux = fti.padding(f,fd,plotQ=False)
	pl.plot(GauXaxis,GauTotal,color='k')
	
	pl.xlabel('Frequency (shifted to %s km/s) (MHz)' % velShift)
	pl.ylabel('%s' % yunits)
	pl.show()


			
			


