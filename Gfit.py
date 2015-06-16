import numpy as np, pylab as pl, G2criteria as G2
from math import *
from scipy.optimize import leastsq
c=2.9979*10**5		#Speed of light, km/s
crlim=0.3
 

def maxabs(arr):
    return arr[np.argmax(np.abs(arr))] 


gauss_fit1 = lambda p, x: p[0]*np.exp(-(x-p[1])**2/(2*(p[2]/2.3548)**2))		#1d Gaussian func
gauss_fit2 = lambda p, x: p[0]*np.exp(-(x-p[1])**2/(2*(p[2]/2.3548)**2))+p[3]*np.exp(-(x-p[4])**2/(2*(p[5]/2.3548)**2)) #2d Gaussian func
e_gauss_fit1 = lambda p, x, y: (gauss_fit1(p,x) -y) #1d Gaussian fit
e_gauss_fit2 = lambda p, x, y: (gauss_fit2(p,x) -y) #2d Gaussian fit


def printFunc(quiet,text,text2='dum',text3='dum',text4='dum',text5='dum',text6 ='dum'):
	if quiet == False:
		if text2!='dum':
			if text3!='dum':
				if text4!='dum':
					if text5!='dum':
						if text6!='dum':
							print text,text2,text3,text4,text5,text6
						else:	print text,text2,text3,text4,text5
					else: print text, text2,text3,text4
				else: print text, text2,text3
			else: print text, text2
		else:  print text

def Best1G(v_a,FMa,v_b,FMb,tBroadF,tNarrowF,xarmin,xarmax,q=False):
	printFunc(q,"In Best1G in Gfit.py, determining which of two Gaussian fits is best.  [H, CR, Wid]1 = ",v_a,"rms1 = %s; [H, CR, Wid]2 = " % (FMa),v_b,"rms1 = %s" % (FMb))
	if list(v_a) == list(v_b):
		printFunc(q, "Two identical fits entered Best1G")
		if FMa != FMb:
			printFunc(q,  "...but somehow they don't have the same values input for FitMeasure.")
		return v_a, FMa
	if v_a[2]> tNarrowF and v_a[2]<tBroadF and abs(v_a[0])<1.E8 and v_a[1]<xarmax and v_a[1]>xarmin:
		if v_b[2]> tNarrowF and v_b[2]<tBroadF and abs(v_b[0])<1.E8 and v_b[1]<xarmax and v_b[1]>xarmin:
			## Both compared fits are "good".  Let's select the one with the lower fit measure.  
			if FMa < FMb:
				return v_a, FMa
			else:
				return v_b, FMb
		else:  
			return v_a, FMa
	else:  
		if v_b[2]> tNarrowF and v_b[2]<tBroadF and abs(v_b[0])<1.E8 and v_b[1]<xarmax and v_b[1]>xarmin:
			return v_b, FMb
		else:  
			## Both fits are "bad"...
			if v_a[1]<xarmin or v_a[1]> xarmax or v_b[1]<xarmin or v_b[1]> xarmax:
				if v_a[1]<xarmin or v_a[1]> xarmax:
					if v_b[1]>xarmax and v_b[1]< xarmax:
						return v_b, FMb
				if v_b[1]<xarmin or v_b[1]> xarmax:
					if v_a[1]>xarmin and v_a[1]< xarmax:
						return v_a, FMa
			if abs(v_a[0])> 1.E8 or abs(v_b[0]) > 1.E8:
				if abs(v_b[0]) < 1.E8:
					return v_b, FMb
				if abs(v_a[0]) < 1.E8:
					return v_a, FMa
				else:   
					printFunc(q,  "both of these have heights that are outrageous.")
					return v_a, FMa
			if v_a[2]< tNarrowF or v_b[2]< tNarrowF:
				if v_b[2] < tNarrowF:
					if v_a[2]>v_b[2]:
						return v_a, FMa
					else:
						return v_b, FMb
				else:  
					return v_b, FMb

			if v_a[2]> tBroadF or v_b[2]> tBroadF:
				if v_b[2] > tBroadF:
					if v_a[2]<v_b[2]:
						return v_a, FMa
					else:
						return v_b, FMb
				else:  
					return v_b, FMb
	printFunc(q, "In Best1G in Gfit.py, determining which of two Gaussian fits is best. Somehow I got to the bottom of Best1G.  Returning v_a, FMa")
	return v_a, FMa
def highlowmid(vv,nn=1):
	
	if vv[1]<vv[4]:
		lowball = vv[1] - nn*abs(vv[2])
		highball = vv[4] + nn*abs(vv[5])
	else: 
		lowball = vv[4] - nn*abs(vv[5])
		highball = vv[1] + nn*abs(vv[2])
	m = (lowball+highball)/2.
	xr=abs(highball-m)
	return xr,m

def plotFunc(plotQ,title,xar,yar,xar1,yar1,GG,col):
	if plotQ == True:
		pl.figure(figsize=(5,3))
		pl.plot(xar, yar,color=col,linewidth=4)
		pl.plot(np.linspace(xar.min(),xar.max(),1000), GG,color='k',linewidth=2)	
		pl.title(title)
		try: pl.axvspan(xar.min(), xar1.min(), color='grey', alpha=0.5)
		except: print "xar1 is empty.  Cannot plot in plotFunc in Gfit.py."
		try: pl.axvspan(xar1.max(), xar.max(), color='grey', alpha=0.5)
		except: print "xar1 is empty.  Cannot plot in plotFunc in Gfit.py."
		pl.xlim([xar.min(),xar.max()])
		pl.show()

def FM(vAR,X,Y,G):
	if len(vAR) == 3:  compMask=abs(X-vAR[1])<0.8*abs(vAR[2])
	if len(vAR) == 6:
		xar_range,mid= highlowmid(vAR,nn=0.8)
		compMask=abs(X-mid)<xar_range
	compY, compG =Y[compMask],G[compMask]
	return np.sqrt(((compY-compG)**2).mean())


def GauFit(xar,yar,rrms,nGau=1,v0=[0.,0.,0.],plotQ=False,tooNarrowv = 2.2, Broadv = 40, tooBroadv=70, tooBroadvABS = 9999,quiet = False, units = 'MHz',compSpread = 20, compWid= 15,DynRange = 100,firstImFreq = 30750.,imFreqRange = 1850.):		##v[0] = height, v[1] = center, v[2] = width (MHz)\
#############################   BinCode = 2 and Bincode = 3 fits should be handled with caution.  Here's the description of these:
#################   binCode = 0:  A 1-component fit was sufficient.
#################   binCode = 1:  A 2-component fit was sufficient.
#################   binCode = 2:  A 2-component fit was deemed necessary, but I couldnt find a good 2-component fit.  I used a 1-component fit. Handle with caution or re-fit manually.
#################   binCode = 3:  The rms in the region fit is significantly higher than that of the spectrum it's coming from.  = bad fit... Handle with caution or re-fit manually.
#################   binCode = 5:  Is WAYY too broad or has an unreasonable line height.  These are in there as placeholders to enable the code to procede (i.e. not continue to loop through finding and attempting to fit the same thing).  You will NEED to refit anything with binCode = 5

	if units == 'GHz':
		xar = xar * 1000
	printFunc(quiet, "  ")
	printFunc(quiet, "  ")
	binCode=0	
	cc='m'
	goTo2=False
	########################### FILL IN THE INITIAL GUESSES FOR THE GAUSSIAN FIT if THEY AREN'T ALREADY INPUT.   ###
	if max(v0)==0:			
		v0[0]=float(maxabs(yar))
		v0[1]=float(xar[np.argmax(np.abs(yar))])
		compWidF = compWid * v0[1]/c
		v0[2]=compWidF
	else: 
		compWidF = compWid * v0[1]/c
	########################### FIND THE IMAGE RMS FROM INPUTS INTO THE FUNCTION.   ###
	try: rsrms = rrms[int((v0[1]-(firstImFreq-imFreqRange/2.))/imFreqRange)]
	except:  
		try: 	rsrms = rrms[int((v0[1]-(firstImFreq-imFreqRange/2.))/imFreqRange)-1]
		except:	rsrms=rrms
	
	########################### DETERMINE THE REASONABLE WIDTHS IN FREQUENCY GIVEN THE INPUTS TO THE FUNCTION.   ###
	if Broadv == 'dum':
		Broadv=40
	if tooBroadv =='dum':
		tooBroadv=70
	if tooBroadvABS=='dum' or tooBroadvABS == 9999:
		tooBroadvABS = tooBroadv
	if tooNarrowv=='dum':
		tooNarrowv=2.2
	BroadF = Broadv * v0[1]/c
	if v0[0]< 0 and tooBroadvABS!=9999:
		tooBroadF = tooBroadvABS * v0[1]/c
		if tooBroadF < BroadF:
			BroadF = .7 * tooBroadF
	else:  
		tooBroadF = tooBroadv * v0[1]/c
	tooNarrowF = tooNarrowv * v0[1]/c
	compSpreadF = compSpread * v0[1]/c
	take1Mask = abs(xar-v0[1])/tooBroadF < .8
	xar1,yar1 = xar[take1Mask],yar[take1Mask]
	printFunc(quiet,"BroadF =  %.3f, tooBroadF = %.3f, tooNarrowF = %.3f, xar1.min() = %.3f, xar1.max() = %.3f" % (BroadF,tooBroadF,tooNarrowF,xar1.min(),xar1.max()))
	########################### FIT A GAUSSIAN GIVEN THE INITIAL GUESSES.   ###
	out = leastsq(e_gauss_fit1, v0[:], args=(xar1, yar1), maxfev=100000, full_output=1)
	v=out[0]
	v[2]=abs(v[2])
	v11G=v
	GG11,GG2 = gauss_fit1(v,np.linspace(xar.min(),xar.max(),1000)), gauss_fit1(v,xar) 
	FitMeasure = FM(v,xar,yar,GG2)
	FM11G=FitMeasure
	printFunc(quiet,'1G initial guesses: [h,cr,wid] = [%.3f,%.3f,%.3f]' % (v0[0],v0[1],v0[2]),'1G 1st fit: [h,cr,wid] = [%.3f,%.3f,%.3f]' % (v[0],v[1],v[2]) )
	printFunc(quiet, "In GauFit, FM11G = ",FM11G)
	plotFunc(plotQ,'1 Gaussian, 1st fit',xar,yar,xar1,yar1,GG11,'b')

	########################### SELECT A REGION TO RE-FIT 1 GAUSSIAN TO.  It will exclude line wings, baseline errors, and blends outside of that range. ###
	printFunc(quiet,"v[2]/tooBroadF = ",v[2]/tooBroadF)
	if v[2]/tooBroadF<1. and v[2]>tooNarrowF:
		RedoMask=abs(xar-v[1])/v[2]<.8	
	else: 
		if v[1] <xar1.max() and v[1]> xar1.min():
			printFunc(quiet, "v[2]/tooNarrowF", v[2]/tooNarrowF)
			printFunc(quiet, "v[0]/rsrms",v[0]/rsrms)
			printFunc(quiet, "v[2]/tooBroadF",v[2]/tooBroadF)
			if v[2]/tooNarrowF>1. and v[0]/rsrms > 10 and v[2]/tooBroadF < 2:
				RedoMask=abs(xar-v[1])/v[2]<.8	
			else:  
				RedoMask=abs(xar-v[1])/tooBroadF<.6		##If the line is broader than what's reasonable, this selects the region closest to the center of the best fit.
		else:
			RedoMask=abs(xar-v0[1])/tooBroadF<.6	##If it got an answer that doesnt really make sense, let's try a more narrow region, closest to the absmax.
	xar1,yar1 = xar[RedoMask],yar[RedoMask]
	if len(xar1) < 5:  
		RedoMask=abs(xar-v[1])/tooBroadF<.6	
		xar1,yar1 = xar[RedoMask],yar[RedoMask]
	try:  printFunc(quiet,"BroadF =  %.3f, tooBroadF = %.3f, xar1.min() = %.3f, xar1.max() = %.3f" % (BroadF,tooBroadF,xar1.min(),xar1.max()))
	except:  
		printFunc(quiet, 'There is still some problem before going in to refit 1Gauss.')
	########################### REFIT   ###
	try:	out = leastsq(e_gauss_fit1, v0[:], args=(xar1, yar1), maxfev=100000, full_output=1)	
	except:  printFunc(quiet, 'The second try of 1G failed.')
	v=out[0] 
	v[2]=abs(v[2])
	GG,GG2 = gauss_fit1(v,np.linspace(xar.min(),xar.max(),1000)), gauss_fit1(v,xar)
	printFunc(quiet, '1G 2nd fit: [h,cr,wid] = [%.3f,%.3f,%.3f]' % (v[0],v[1],v[2])) 
	plotFunc(plotQ,'1 Gaussian, 2nd fit',xar,yar,xar1,yar1,GG,'b')
	v1G=v			
	FitMeasure = FM(v,xar,yar,GG2)
	FM1G= float(FitMeasure)
	########################## IF THE SECOND 1G FIT LOOKS REASONABLE, TAKE THE BETTER OF THE 1ST & 2ND FITS.   ###
	if abs(v[0])<1.E8 and abs(v[2])<tooBroadF and abs(v[2])>tooNarrowF and v[1]<xar1.max() and v[1]>xar1.min():	##Then 1GAUSS looks reasonable.  I WOULD LIKE TO CHANGE v[2]<20 TO 12.   ###
		printFunc(quiet, "Got into the 1Gauss isn't totally egregious")
		vBEST1,FMBEST1 = Best1G(v11G,FM11G,v1G,FM1G,tooBroadF,tooNarrowF,xar.min(),xar.max(),q=quiet)
		v, FMBEST1 = Best1G(v11G,FM11G,v1G,FM1G,tooBroadF,tooNarrowF,xar.min(),xar.max(),q=quiet)
	########################## IF THE SECOND 1G FIT DOESNT LOOK REASONABLE, TRY FITTING A DIFFERENT REGION (A BROADER REGION) AND THEN TAKE THE BEST OF THE 1ST, 2ND, AND NEW FIT   ###
	else: 
		printFunc(quiet, "The second iteration of 1Gauss isn't convincing.")
		take3Mask = abs(xar-v11G[1])/tooBroadF < 1  ## Try just fitting to a bigger region".
		xar1,yar1 = xar[take3Mask],yar[take3Mask]
		if len(xar1) == 0:  
			take3Mask = abs(xar-v0[1])/tooBroadF < 1  ## Try just fitting to a bigger region".
			xar1,yar1 = xar[take3Mask],yar[take3Mask]
		try:  out = leastsq(e_gauss_fit1, v0[:], args=(xar1, yar1), maxfev=100000, full_output=1)		##This fits the region set by the 1st 1G fit. 
		except: 
			printFunc(True, 'The 2nd 1 G was bad, the 3rd try didnt work.')
		v = out[0]
		v[2]=abs(v[2])
		printFunc(quiet, '1G 2nd fit was unreasonable, new fit: [h,cr,wid] = [%.3f,%.3f,%.3f]' % (v[0],v[1],v[2])) 
		plotFunc(plotQ,'1 Gaussian, new fit',xar,yar,xar1,yar1,GG,'b')
		FitMeasure = float(FM(v,xar,yar,GG2))
		## Find the better fit between the 3rd try and the second.  
		vbetter,FMbetter = Best1G(v,FitMeasure,v1G,FM1G,tooBroadF,tooNarrowF,xar.min(),xar.max(),q=quiet)
		## Find the better fit between the winner of round 1 and the second.
		vBEST1,FMBEST1 = Best1G(vbetter,FMbetter,v11G,FM11G,tooBroadF,tooNarrowF,xar.min(),xar.max(),q=quiet)
		v, FMBEST1 = Best1G(v11G,FM11G,v1G,FM1G,tooBroadF,tooNarrowF,xar.min(),xar.max(),q=quiet)
		if v[2]< tooNarrowF or v[2]>tooBroadF or abs(v[0])>1.E8 or v[1]>xar.max() or v[1]<xar.min():
			printFunc(quiet, "The most reasonable 1G fit is still unreasonable.  binCode = 3.")
			binCode = 3
			cc='g'
	########################### DETERMINE IF 2GAUSS IS CALLED FOR OR 1GAUSS IS SUFFICIENT ############################	
	printFunc(quiet," Passing these parameters to G2.criteria:  cr-cr0 = %.3f, FM = %.3f, ImageRMS =  %.3f, height = %.3f, absMax is at %.3f" % (vBEST1[1]-v0[1],FMBEST1, rsrms, vBEST1[0],v0[0]))
	if G2.criteria(v[1],v0[1],v[2],FMBEST1,v[0],v0[0],rsrms,BroadF, q=quiet)==True:	
			goTo2='True'
			binCode=1
			cc = 'y'
			printFunc(quiet, "MET the 2Gauss criterion the first try!")
	
	else: printFunc(quiet, "didnt meet the 2Gauss criterion the first try!")	
	SNR = abs(v[0]/rsrms)
	FMrms= FitMeasure/rsrms
	
	printFunc(quiet, "!!!!!!Best 1 Gaussian fit obtained for Freq =", v0[1])
	printFunc(quiet, "!!!       SNR/DynRange * 1/FMrms = ",SNR/DynRange * 1/FMrms)
	printFunc(quiet, "!!! FM/rms = ",FMrms )
	printFunc(quiet, "!!! SNR/DynRange = ",SNR/DynRange )
##########################################################################################
####### DONE IF ONLY 1GAUSS IS NEEDED.  IF 2GAUSS, NOW MOVING INTO FITTING 2GAUSS. #######
##########################################################################################
	if goTo2=='True':	
		printFunc(quiet, "got into 2Gaus ")
		printFunc(quiet, "v start = ",v[0],v[1],v[2])
		nGau = 2

############################### FILL THE INPUT GUESS ARRAY & ESTABLISH A REGION TO FIT THAT CONTAINS IT.  IT ACTUALLY TRIES TO USE THE REGION FROM ABOVE, AND ONLY CHANGES xar1 IF IT NEEDS TO.  THEN FIT A 2GAUSS PROFILE.  ###
		G2_input=[v0[0],v0[1],compWidF]	
		printFunc(quiet, "Best 1Gauss is centered at = ", vBEST1[1])
		direction = (v[1]-v0[1])/abs(v[1]-v0[1])
		FrangeRequired = BroadF/10

         ########################### PUT IN SOME INITIAL GUESSES.   ###
		if v0[1]+ direction * compSpreadF > xar.min() and v0[1]+ direction * compSpreadF < xar.max():
		  # we expect to be able to use our best first guesses.
			G2_input.extend([G2_input[0]*0.6,G2_input[1]+2.3*(v[1]-v0[1]),compWidF])
		else:  	
			if G2_input[1]+2.3*(vBEST1[1]-v0[1])> xar.min() and G2_input[1]+2.3*(vBEST1[1]-v0[1])< xar.max():
				G2_input.extend([G2_input[0]*0.6,G2_input[1]+2.3*(v[1]-v0[1]),compWidF])
			else:  
				printFunc(quiet, "got in IF place in fill input guess array", "G2_input[1]-1.*(v[1]-v0[1]) =",G2_input[1]-1.*(v[1]-v0[1])) 
				printFunc(quiet, "G2_input[1]", G2_input[1], "v[1]",v[1],"v0[1]",v0[1])
				direction = -1*direction
				G2_input.extend([G2_input[0]*0.6,G2_input[1]+ direction * compSpreadF,compWidF])


         ########################### ADJUST THE FIT INPUT AND THE REGION WE'RE FITTING.   ###
		if abs(G2_input[4]-G2_input[1])/compSpreadF<.6 or abs(G2_input[4]-G2_input[1])/compSpreadF>2.5:	
		## Then our initial guesses aren't reasonable given the typical spread in components.
			G2_input[4]=G2_input[1]+compSpreadF * direction
		if G2_input[4]-xar1.min() < FrangeRequired  or xar1.max() - G2_input[4] < FrangeRequired: 
			printFunc(quiet, "Re-setting a region")
			if G2_input[1]- xar.min() < FrangeRequired or xar.max() - G2_input[1] < FrangeRequired:
				print "Something's super wrong in GauFit!"
			if G2_input[4] - xar.min() < FrangeRequired or xar.max() - G2_input[4] < FrangeRequired:
				if direction == 1:
					xlim = xar.max()
				if direction == -1:
					xlim = xar.min()
				G2_input[4]= xlim - direction * FrangeRequired
			mFin = np.mean([G2_input[1],G2_input[4]])
			if 1.3*abs(G2_input[1]-mFin) / v[2] < 1.5:
				doMask=abs(xar-np.mean([G2_input[1],G2_input[4]]))/abs(G2_input[2])<1.5		
			else:  
				doMask=abs(xar-np.mean([G2_input[1],G2_input[4]]))/abs(G2_input[1]-mFin)<1.3
			xar1,yar1=xar[doMask],yar[doMask]	

		if len(xar1) < 5:   ##Somehow... I think I've covered all my bases, but just in case		
			xar1,yar1 = xar[take1Mask],yar[take1Mask]
			if len(xar1) < 5:
				binCode = 2
				goTo2 = 'getOut'  ## Give it a way out.
				printFunc(quiet, "************************")
				printFunc(quiet, " ")
				printFunc(quiet, "Somehow, I can't get a good combination for the input guesses and regions for 2Gauss.  Am continuing with 1 Gauss.")
				printFunc(quiet, " ")
				printFunc(quiet, "************************")
			
		printFunc(quiet, '2G initial guesses: G1 [h,cr,wid], G2 [h,cr,wid] = [%.3f, %.3f, %.3f], [%.3f,%.3f,%.3f]' % (G2_input[0],G2_input[1],G2_input[2],G2_input[3],G2_input[4],G2_input[5]))
		printFunc(quiet, 'Fitting data using the frequency range [%.3f, %.3f]' % (xar1.min(),xar1.max()))
         ########################### DONE FILLING THE INPUT GUESS ARRAY & SELECTING THE REGION. ###

         ########################### FIT THE 1ST TAKE ON THE 2GAUSS FIT.   ### 
		try:  out = leastsq(e_gauss_fit2, G2_input[:], args=(xar1, yar1), maxfev=100000, full_output=1) 	#Gauss Fit
		except:  
			goTo2='getOut'
			binCode=2
			cc = 'g'
			nGau = 1
			printFunc(quiet, "2Gauss didnt go so well")
			v=vBEST1
			FitMeasure = FMBEST1
	if goTo2=='True':
		v=out[0] #fit parameters out
		v[2] = abs(v[2])
		v[5] = abs(v[5])
		printFunc(quiet, 'first try, 2G: G1 [h,cr,wid], G2 [h,cr,wid] = [%.3f,%.3f,%.3f], [%.3f,%.3f,%.3f]' % (v[0],v[1],v[2],v[3],v[4],v[5]))
		GG,GG2=gauss_fit2(v,np.linspace(xar.min(),xar.max(),1000)), gauss_fit2(v,xar)
		plotFunc(plotQ,'First 2Gauss fit.',xar,yar,xar1,yar1,GG,cc)

############################### GENERATE A SECOND 2GAUSS FIT.  SELECT A REGION TO RE-FIT AND TRY FITTING, WITH THE SAME INPUT PARAMETERS FROM ABOVE, BUT WITH A DIFFERENT REGION. ###

         ########################### ADJUST THE REGION WE'RE FITTING.   ###
		xar_range,mid= highlowmid(v,nn=0.8)
		printFunc(quiet, "xar_range,mid =",xar_range,mid)
		if abs(v[0])<1.E8 and abs(v[3])<1.E8 and v[1]<xar.max() and v[1]>xar.min() and v[4]<xar.max() and v[4]>xar.min(): 
			if xar_range > tooNarrowF and xar_range< tooBroadF*1.5 and v[2]<tooBroadF and v[5]<tooBroadF:  
				doOverMask=abs(xar-mid)<xar_range
				printFunc(quiet, "xar_range = %0.3f" % xar_range)
				printFunc(quiet, "set doOverMask in else, if")
			else:
				GinpRange,m = highlowmid(G2_input,0.7)
				printFunc(quiet, "GinpRange,m =",GinpRange,m)
				doOverMask=abs(xar-m)<GinpRange
				printFunc(quiet, "set doOverMask in else, else")
		else: 
			printFunc(quiet, "NEED TO ADJUST THE REGION SIGNIFICANTLY")			    
			if abs(v[2])<tooBroadF or abs(v[5])<tooBroadF:		##Might be able to salvage things using one of the 2-gauss components if one has a decent parameter for width.  Also, we'll populate an input array for a 1component fit if we can't find a good 2comp one.
				if abs(v[2])<abs(v[5]):
				
				    if v[1] > xar.min() and v[1] < xar.max() and v[0]/v0[0] > 0 and abs(v[0])< 1.E8:
					v[2]=abs(v[2])
					doOverMask = abs(xar-v[1])<.7*v[2]	##Make it a slightly more narrow range.
				    else:  	
					if v[4] > xar.min() and v[4] < xar.max():
					    v[5]=abs(v[5])
					    doOverMask = abs(xar-v[4])<.7*v[5]
					else:  
						printFunc(quiet, " DIDN'T FIND A DECENT ADJUSTMENT>>>")
						goTo2='getOut'
						binCode=2
						cc = 'g'
						nGau = 1
						printFunc(quiet, "2Gauss didnt go so well.  Will use a 1G fit.")
						v=vBEST1
						FitMeasure = FMBEST1
				else:
				    if abs(v[5])<abs(v[2]): 
				    	if v[4] > xar.min() and v[4] < xar.max() and v[3]/v0[0] > 0 and abs(v[3])< 1.E8: 
				 	    v[5]=abs(v[5])
					    doOverMask = abs(xar-v[4])<.7*v[5]
					else:  	
					    if v[1] > xar.min() and v[1] < xar.max():
					        v[2]=abs(v[2])
					        doOverMask = abs(xar-v[1])<.7*v[2]	
					    else:  
						printFunc(quiet, " DIDN'T FIND A DECENT ADJUSTMENT>>>")
						goTo2='getOut'
						binCode=2
						cc = 'g'
						nGau = 1
						printFunc(quiet, "2Gauss didnt go so well.  Will use a 1G fit.")
						v=vBEST1
						FitMeasure = FMBEST1
			else:  
				printFunc(quiet, " DIDN'T FIND A DECENT ADJUSTMENT>>>")
				goTo2='getOut'
				binCode=2
				cc = 'g'
				nGau = 1
				printFunc(quiet, "2Gauss didnt go so well.  Will use a 1G fit.")
				v=vBEST1
				FitMeasure = FMBEST1

	if goTo2=='True':
		xar1,yar1 =xar[doOverMask],yar[doOverMask]
		printFunc(quiet, "min, max of the region we're re-fitting in 2Gauss = %.2f, %.2f" % (xar1.min(),xar1.max()))


         ########################### FIT THE 2ND TAKE ON THE 2GAUSS FIT.   ### 
		try:
			out = leastsq(e_gauss_fit2, G2_input[:], args=(xar1, yar1), maxfev=100000, full_output=1) 	#Gauss Fit
		except:
			printFunc(quiet, "the second Gaussian fit didnt work")
		v=out[0] #fit parameters out
		printFunc(quiet, 'best fit, 2G: G1 [h,cr,wid], G2 [h,cr,wid] = [%.3f,%.3f,%.3f], [%.3f,%.3f,%.3f]' % (v[0],v[1],v[2],v[3],v[4],v[5]))
		v[2]=abs(v[2])
		v[5]=abs(v[5])
		### DONE FITTING 2 GAUSS ### 

		### EVALUATE THE RMS RESIDUAL ###
		GG, GG2=gauss_fit2(v,np.linspace(xar.min(),xar.max(),1000)), gauss_fit2(v,xar) 
		xar_range,mid= highlowmid(v,nn=0.8)
		FitMeasure = FM(v,xar,yar,GG2)
		plotFunc(plotQ,'2nd 2Gauss fit.',xar,yar,xar1,yar1,GG,cc)
		### DONE EVALUATING THE RMS RESIDUAL ###

############################### BUILD A 'SOMETHING DOESNT WORK PLACE' IN WHICH I CAN HAVE ONE MORE GO AT PRODUCING A DECENT 2GAUSS FIT OR DECIDE TO USE THE BEST 1G FIT.   ###

         ###########################  TRY A 2-GAUSS FIT ON A BROADER AREA.  IF THAT DOESN'T WORK, THEN TAKE THE 1ST 2GAUSS FIT IF IT WAS OKAY.  IF IT WASN'T OKAY, TAKE THE BEST 1GAUSS FIT.   ###
		FMrms= FitMeasure/rsrms
		printFunc(quiet, "After 2nd try at 2Gauss, FitMeasure/rsrms= % 0.2f" % FMrms) 
		if v[0] > v[3]:
			SNR = abs(v[0]/rsrms)
		else:  
			SNR = abs(v[3]/rsrms)
		printFunc(quiet,  "!!!!!!!!! SNR/DynRange * 1/FMrms = ", SNR/DynRange * 1/FMrms)
		printFunc(quiet,  "!!!!      FM/rms = ",FMrms)
		if v[2]/tooNarrowF<.6 or v[5]/tooNarrowF<0.6 or abs(v[0])>1.E8 or abs(v[3])>1.E8 or v[2]/tooBroadF>1 or v[5]/tooBroadF>1 or v[1]>xar1.max() or v[1]<xar1.min() or v[4]>xar1.max() or v[4]<xar1.min(): 
			printFunc(quiet, "I got into the 2Gauss fit -- SOMETHING DOESNT WORK PLACE.")
			plotFunc(plotQ,'Failed 2Gauss fit.',xar,yar,xar1,yar1,GG,cc)

			if v11G[2]<tooBroadF*1.3 and v11G[2]>tooNarrowF and v11G[1] > xar1.min() and v11G[1] > xar1.max():
				RedoMask=abs(xar-v11G[1])/v11G[2]<.8		
			else: 
				if v11G[1] > xar1.min() and v11G[1] > xar1.max():
					RedoMask=abs(xar-v11G[1])<BroadF
				else: 
					RedoMask=abs(xar-v0[1])<BroadF
			xar1,yar1=xar[RedoMask],yar[RedoMask]
			### DONE SELECTING THE REGION TO FIT ### 
	
			### NOW FIT AND EVALUATE FIT.   ### 
			try:  out = leastsq(e_gauss_fit2, G2_input[:], args=(xar1, yar1), maxfev=100000, full_output=1) 	#Gauss Fit 
			except:  
				goTo2='getOut'
				binCode=2
				cc = 'g'
				nGau = 1
				printFunc(quiet, "2Gauss re-working FAILED")
				v=vBEST1
				FitMeasure = FMBEST1
			if goTo2!='getOut':			
			   v=out[0] #fit parameters out
			   v[2] = abs(v[2])
			   v[5] = abs(v[5])
			   printFunc(quiet, 'After 2G didnt work: G1 [h,cr,wid], G2 [h,cr,wid] = [%.3f,%.3f,%.3f], [%.3f,%.3f,%.3f]' % (v[0],v[1],v[2],v[3],v[4],v[5]))
			   GG, GG2=gauss_fit2(v,np.linspace(xar.min(),xar.max(),1000)), gauss_fit2(v,xar)
			   FitMeasure = FM(v,xar,yar,GG2)
			   plotFunc(plotQ,'2Gauss re-fit after first try didnt work.',xar,yar,xar1,yar1,GG,cc)
			   if v[2]<tooNarrowF or v[5]<tooNarrowF or abs(v[0])>1.E8 or abs(v[3])>1.E8 or v[2]/tooBroadF>1 or v[5]/tooBroadF>1 or v[1]>xar1.max() or v[1]<xar1.min() or v[4]>xar1.max() or v[4]<xar1.min(): 
				FitMeasure=9999
				nGau = 1
				binCode=2
				cc = 'g'
				if abs(v[2])<tooBroadF or abs(v[5])<tooBroadF:		##Might be able to salvage things using one of the 2-gauss components if one has a decent parameter for width.
					if abs(v[2])<abs(v[5]):
					    if v[1] > xar.min() and v[1] < xar.max() and v[3]/v0[0] > 0:
						v[2]=abs(v[2])
						doOverMask = abs(xar-v[1])<.7*v[2]	##Make it a slightly more narrow range.
						vINP = v[0:3]	
						vINP[0]= v0[0]
					    else:  	
						if v[4] > xar.min() and v[4] < xar.max():
							v[5]=abs(v[5])
							doOverMask = abs(xar-v[4])<.7*v[5]
							vINP = v[3:6]	
							vINP[0]= v0[0]
						else:  
							v[5]=abs(v[5])
							doOverMask=abs(xar-v[4])<Broadv*.4
							vINP = v[3:6]	
							vINP[0]= v0[0]

					if abs(v[5])<abs(v[2]): 
					    if v[4] > xar.min() and v[4] < xar.max() and v[3]/v0[0] > 0: 
						v[5]=abs(v[5])
						doOverMask = abs(xar-v[3])<.7*v[5]
						vINP = v[3:6]
						vINP[0]= v0[0]	
					    else:  
						if v[1] > xar.min() and v[1] < xar.max():	
							v[2]=abs(v[2])
							doOverMask = abs(xar-v[1])<.7*v[2]
							vINP = v[0:3]
							vINP[0]= v0[0]	
						else:  
							v[2]=abs(v[2])
							doOverMask=abs(xar-v[1])<Broadv*.4
							vINP = v[0:3]	
							vINP[0]= v0[0]	

					xar1, yar1 =xar[doOverMask],yar[doOverMask] 
					printFunc(quiet,  "vINP = [%.3f, %.3f, %.3f]" %(vINP[0],vINP[1],vINP[2]))
					try:  out = leastsq(e_gauss_fit1, vINP[:], args=(xar1, yar1), maxfev=100000, full_output=1)
					except:  FitMeasure = 9999
					v=out[0]
					GG, GG2=gauss_fit1(v,np.linspace(xar.min(),xar.max(),1000)), gauss_fit1(v,xar)
					FitMeasure = FM(v,xar,yar,GG2)
					plotFunc(plotQ,'1Gauss re-fit after being unable to get a good 2G fit.',xar,yar,xar1,yar1,GG,cc)
					printFunc(quiet, "1Gauss re-fit, v = ",v)
					vBEST1,FMBEST1 = Best1G(vBEST1,FMBEST1,v,FitMeasure,tooBroadF,tooNarrowF,xar.min(),xar.max(),q=quiet)
					v, FitMeasure = Best1G(vBEST1,FMBEST1,v,FitMeasure,tooBroadF,tooNarrowF,xar.min(),xar.max(),q=quiet)
					
	if FitMeasure/rsrms > 2 and abs(FitMeasure/v[0]) > .1:  ###This is a bad fit.  
		if len(v)==6:
		    if abs(FitMeasure/v[3]) >.1:
			binCode= 3
			cc = 'g'
		if len(v)==3:
			binCode= 3
			cc = 'g'
	if v[2]/tooBroadF > 5 or tooNarrowF/v[2] > 5 or abs(v[0])/abs(v0[0])>2:  ###This is a bad fit.  
			binCode= 5
			cc = 'r'
			v=v0
			nGau = 1
			FitMeasure = 9999
	if binCode==1:
		if abs(v[3])/abs(v0[0])>2:  ###This is a bad fit.  
			binCode= 5
			cc = 'r'
			v=v0
			nGau = 1
			FitMeasure = 9999
			
	if plotQ==True:
		try: GG=gauss_fit2(v,np.linspace(xar.min(),xar.max(),1000))
		except: GG=gauss_fit1(v,np.linspace(xar.min(),xar.max(),1000))
		plotFunc(plotQ,'Final Fit',xar,yar,xar1,yar1,GG,cc)
	if nGau == 1:
		printFunc(quiet, "height = %.3f, cr = %.3f, width = %.3f, FitMeasure = %.3f, binCode = %i " % (v[0],v[1],v[2],FitMeasure,binCode))
		return [binCode,v[0],v[1],v[2],FitMeasure,v[1]-v0[1]]
	if nGau == 2:
		printFunc(quiet, "height = %.3f, cr = %.3f, width = %.3f, height = %.3f, cr = %.3f, width = %.3f, FitMeasure = %.3f, binCode = %i" % (v[0],v[1],v[2],v[3],v[4],v[5],FitMeasure,binCode))
		return [binCode,v[0],v[1],v[2],v[3],v[4],v[5],FitMeasure]

