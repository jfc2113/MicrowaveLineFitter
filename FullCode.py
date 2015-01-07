import numpy as np, pylab as pl, LineFind as LF, Gfit_dynRange_4 as Gfit, LineID as lID, FTinterp as fti, PBcsv
from math import *
c=2.9979E5

thresh = 3.5
wayout='False'  ## If wayout == 'True', after determining the regions with lines, the code will plot the full spectrum and the line regions without fitting anything.  
annotateQ='True'
Compcsv='False'
Outcsv = 'False'
GoutPlotQ = False
GQuiet = True

#Compcsv='True'
#Outcsv = True
#wayout = True
#GoutPlotQ=True
GQuiet = False


#regions=['LMH','h']
#regions =['h_dum']
#regions=['N','LMH','L']
regions=['L','K4']#,'K4','L']#,'K4']

regions=['K6']#,'K4','L']#,'K4']
compFile1 = 'recomb.csv'
compFile2 = 'DetectedLines2.csv'
def OF(r):
	return "ALLFITS_%s_Sept14_4.csv" % r

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
		#if Goutp[0]!=2:
		if Goutp[0] == 2:
			print "GOUTp = ",Goutp
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
		#if Goutp[0]!=2:
		if Goutp[0] == 2:
			print "GOUTp = ",Goutp
		GT+=Gfit.gauss_fit2(Goutp[1:7],bigFreq)
	return lineAr,hAr,crAr,widAr,FMAr,croffAr,binAr,GG,GT

for r in regions:
	outfile = OF(r)
	MoreGs=0
	fig = pl.figure()
	ax = fig.add_subplot(111)
	lineCounts=0
	f,fd = LF.spect('HP_%s_ALL.txt' %(r))
#	f,fd = LF.spect('%s_all.txt' %(r))
	freq,flux,regRMS=LF.LineFind(f,fd,plotQ='False',rmsthresh1=thresh)
	print "wayout!!!", wayout
	if wayout ==True:
	   print "see, wayout = True"
	   pl.plot(f,fd,'c',linewidth=5)
	   starting = 0
	   freqDif=freq[1:]-freq[:len(freq)-1]
	   for iiii in range((len(freqDif))):
		if abs(freqDif[iiii])>1.9:
			ending = iiii+1
#			print starting, ending
			if ending-starting > 3:
				smallFreq=freq[starting:ending]	
				smallFlux=flux[starting:ending]
				pl.plot(smallFreq,smallFlux,color='y',linewidth=5)
				interpFreq,interpFlux=fti.padding(smallFreq,smallFlux,plotQ='False')
				pl.plot(interpFreq,interpFlux,color='b',linewidth=1)
			if abs(smallFreq.mean() - 38271) < 10 or abs(smallFreq.mean() - 33704) < 10 or abs(smallFreq.mean() - 36780) < 40:  
				print smallFreq.mean(),":",  starting, ending
			starting = ending

  	   pl.xlabel('Frequency (MHz)')
	   pl.ylabel('Brightness (mJy/arcsec$^{-2}$)')
	   pl.show()
	   break
	GauTotal = np.zeros((100000))
	GauXaxis = np.linspace(f.min(),f.max(),100000)
	tellmewtfmask=abs(freq-46256)<100
	regRMS0=regRMS
#	if r =='h' or r=='LMH' or r=='N':
#		regRMS=4.
	starting = 0
	freqDif=freq[1:]-freq[:len(freq)-1]
	hs=[]
	FMs=[]
	croffs=[]
	crs=[]
	wids=[]
	binCodes=[]
	lineIDs=[]
	for iiii in range((len(freqDif))):
		if abs(freqDif[iiii])>1.9:
			ending = iiii+1
		#	print "iiii= ", iiii
#			print starting, ending
			if ending-starting > 3:
########  NEXT LINE = UNCOMMENT IF YOU WANT TO SHOW SPECIFIC PARTS OF THE SPECTRUM   ########
		#	   if starting > 2220 and starting<2240: 	
		#	   if starting > 850 and starting<880: 		###In case you want to target specific parts of the spectrum	
		#	   if starting > 700 and starting<740: 	
		#	   if starting > 1120 and starting<1150: 
		#	   if starting > 4400 and starting<5000: 
	#		   if starting > 330 and starting<380: 				
				print "iiii= ", iiii
				smallFreq=freq[starting:ending]	
				smallFlux=flux[starting:ending]
				interpFreq,interpFlux=fti.padding(smallFreq,smallFlux,plotQ='False')
				print " \n \n \n"
				print "In FullCode:  Passing interpFreq.min() = %.2f to interpFreq.max() =%.2f to Gfit" % (interpFreq.min(),interpFreq.max())
				Gout=Gfit.GauFit(interpFreq,interpFlux,regRMS,1,[0,0,0],quiet=GQuiet, plotQ=GoutPlotQ,tooBroadABS=30)	## Broad absorption is ~5MHz wide.
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
					print "Fitting G in iteration %i" % dumi
					interpFreq,interpFlux=fti.padding(smallFreq,smallFlux,plotQ='False')
					Gout=Gfit.GauFit(interpFreq,interpFlux,regRMS,1,[0,0,0],quiet=GQuiet,plotQ=GoutPlotQ,tooBroadABS=35)
					print "Gout[0]",Gout[0]
					lineIDs,hs,crs,wids,FMs,croffs,binCodes,GGG,GauTotal = listfill(lineIDs,hs,crs,wids,FMs,croffs,binCodes,lineCounts,Gout,smallFreq,GauTotal,GauXaxis)
					lineCounts+=1
			starting = ending

	smallFreq=freq[starting:len(freq)]	
	smallFlux=flux[starting:len(freq)]
	interpFreq,interpFlux=fti.padding(smallFreq,smallFlux,plotQ='False')
	print "  "
	print "In FullCode:  Passing interpFreq.min() = %.2f to interpFreq.max() =%.2f to Gfit" % (interpFreq.min(),interpFreq.max())
	Gout=Gfit.GauFit(interpFreq,interpFlux,regRMS,1,[0,0,0],quiet=GQuiet,plotQ=GoutPlotQ,tooBroadABS=35)
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
#		print "Fitting G in iteration %i" % dumi
		interpFreq,interpFlux=fti.padding(smallFreq,smallFlux,plotQ='False')
		Gout=Gfit.GauFit(interpFreq,interpFlux,regRMS,1,[0,0,0],quiet=True,plotQ=GoutPlotQ)
		print "Gout[0]",Gout[0]
		lineIDs,hs,crs,wids,FMs,croffs,binCodes,GGG,GauTotal = listfill(lineIDs,hs,crs,wids,FMs,croffs,binCodes,lineCounts,Gout,smallFreq,GauTotal,GauXaxis)
		lineCounts+=1
		
	regionMaster=zip(lineIDs,crs,croffs,hs,wids,FMs,binCodes)
	crs=np.asarray(crs)
	if Outcsv == True or Outcsv == 'True':
		lID.lineCompare('recomb.csv','recombFits_%s.csv' % r,regionMaster,crs,plotQ='False',writeout=True)
		lID.lineCompare('DetectedLines2.csv','molecFits_%s.csv' % r,regionMaster,crs,plotQ='False',writeout=True)
	#	lID.lineCompareAll('recomb.csv','ALLfits_%s_MayDum.csv' % r,regionMaster,plotQ='False',writeout=True,infile2 = 'DetectedLines2.csv')
		lID.lineCompareAll(compFile1,outfile,regionMaster,plotQ='False',writeout=True,infile2 = compFile2)
		print "should have written outfile = %s out" % outfile


	for icrfreqi in range(len(regRMS)):
		crf = 30750+1850*icrfreqi
		rmsNow= regRMS[icrfreqi]
		rmsPB= PBcsv.PBcor(r,crf,rmsNow)
		print "for image at crFreq = ",crf, "RMS = ",rmsNow , "After PBcor:", rmsPB
	pl.plot(f,fd,color = 'c')#,linewidth=4)
	interpFreq,interpFlux = fti.padding(f,fd,plotQ=False)
	pl.plot(GauXaxis,GauTotal,color='k')
	
	pl.xlabel('Frequency (shifted to 64km/s) (MHz)')
	pl.ylabel('Brightness (mJy/arcsec$^{-2}$)')
	if Compcsv=='True':
		lID.lineCompareAll(compFile1,outfile,regionMaster,plotQ='True',writeout=False,infile2 = compFile2)
	pl.show()


			
			


