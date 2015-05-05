#####################################################################################################
##############################  This is the one you actually run!      ##############################
##############################  See VelEDIT.py for all parameters.     ##############################
#####################################################################################################


#####################################################################################################
#####################  This script will conver the output of FullCode.py to a concise csv file with 
#####################  line IDs and gaussian fits in velocity space.
#####################  Before running, re-fit any lines in ALLFITS with bin = 5, and fix any other 
#####################  fits as desired.  Edit material in the file VelocityEditFile.py.
#####################  Also, you can optionally prepare 1-2 csv files containing prioritized line IDs.  
#####################  As it is currently written, the code inspects assigned recombination lines (priority 1)
#####################  for kinematic consistency first.  The script does not require a csv file input 
#####################  for the recombination lines at this stage however.  Then (priority 2) it looks for 
#####################  lines that you specify as expecting to be strong, in a file "strongMol.csv".
#####################  Third, it handles  foreground absorption by very strong lines held in 'SAClines.csv'
#####################  (SAC stands for spiral arm cloud).  Fourth, it assigns lines detected in ALLFITS.csv.  
#####################  Finally, it handles unidentified transitions.
#####################  Users can specify either 1 or 2 csv files to establish the priorities as described below.


###################################################################################################################

#####################   DEFINE FUNCTIONS AND RUN THROUGH A LOOP USING INFORMATION PROVIDED IN THE VelocityEditFile.py. 
import numpy as np, pylab as pl, csv
import VelocityEditFile as VE
c=2.99792E5		## Speed of light in km/s


import numpy as np, pylab as pl, csv


def lineGrab(infile,units='GHz'):
	try: tt=csv.reader(open(infile,'rt'),delimiter=':')
	except: print "couldn't open the file %s" % infile
	fs,IDs,trans=[],[],[]
	for row in tt:
		try: 	f=float(row[0])
		except:	continue
		if units == 'GHz':	fs.append(f*1000)
		if units == 'MHz':	fs.append(f)
		if units == 'kHz':	fs.append(f/1000.)	
		IDs.append(row[1])
		trans.append(row[2])		
	return np.asarray(fs),  np.asarray(IDs), np.asarray(trans)
	

				
def ArrayFiller(reg,ALLFITSfile,coolbins =[0,1,2,3,4],fUnitAF='GHz',fUnitCSV='GHz'):
	t=csv.reader(open(ALLFITSfile,'rt'),delimiter=':')
	try: coolbins = VE.goodBins(reg)
	except: print "the good bins are assumed to be 0,1,2,3, and 4"
	try: strongFile,SACfile,fUnitCSV = VE.strongAndSAC(reg)
	except: print "Didn't get anything from VelocityEditFile.py"
	try: fSAC,IDSAC,transSAC = lineGrab(SACfile,fUnitCSV)
	except: 
		print "Could not get data from SACfile.  May not have set a file name.  Proceeding."
		fSAC,IDSAC,transSAC = np.asarray([]),np.asarray([]),np.asarray([])
	try: fStrong,IDStrong,transStrong = lineGrab(strongFile,fUnitCSV)
	except: 
		print "Could not get data from StrongLines.  May not have set a file name.  Proceeding."
		fStrong,IDStrong,transStrong = np.asarray([]),np.asarray([]),np.asarray([])
	Umask= IDStrong=='U'		## in case you have known 'U' lines in one of the strong line files
	try: fStrongU=fStrong[Umask]
	except: fStrongU=np.asarray([])
	ffits,flines,linetype,heights,wids,transitions,speciess=[],[],[],[],[],[],[]
	iii = 0
	for row in t:	
		try: f=float(row[1])			# This should be in MHz, by the default output of FullCode.py.
		except: continue
		if int(row[6]) not in coolbins:		# If it was a clearly bad fit
			print "Ignoring a bad fit at freq = %.2f MHz" % f
			continue
		ffits.append(f)
		heights.append(float(row[3]))
		wids.append(float(row[4]))
		try: fline=float(row[9])	#I had it in GHz.  Currently plotting in MHz.
		except: 
			iii+=1
			if len(fStrongU)> 0:
			   if np.min(abs((fStrongU-f)/f*c))<10:		#If the line is within 10km/s of a u-line.
				linetype.append('StrongU')	
				flines.append(fStrongU[np.argmin(abs(fStrongU-f))])
				transitions.append('Known U-line')
				speciess.append('Known U-line')
			   else:
				linetype.append('unidentified')
				flines.append(0.0)
				transitions.append('un')
				speciess.append('un')
			else:
				linetype.append('unidentified')
				flines.append(0.0)
				transitions.append('un')
				speciess.append('un')
			continue
		if fUnitAF=='GHz': 	flines.append(1000*fline)
		if fUnitAF=='MHz': 	flines.append(fline)
		if fUnitAF=='kHz': 	flines.append(fline/1000.)
		
		if row[12]=='Recomb':
			if row[7]=='H&alpha;':
				linetype.append('H_alpha')
				transitions.append(row[10])
				speciess.append(row[7])
			if row[7]=='He&alpha;':
				linetype.append('He_alpha')
				transitions.append(row[10])
				speciess.append(row[7])
			if row[7]=='He&beta;':
				linetype.append('He_beta')
				transitions.append(row[10])
				speciess.append(row[7])
			if row[7]=='H&beta;':
				linetype.append('H_beta')
				transitions.append(row[10])
				speciess.append(row[7])
			if row[7]=='H&gamma;':
				linetype.append('H_gamma')
				transitions.append(row[10])
				speciess.append(row[7])
			if row[7]=='H&delta;':
				linetype.append('H_delta')
				transitions.append(row[10])
				speciess.append(row[7])
			if row[7]=='H&epsilon;':
				linetype.append('H_epsilon')
				transitions.append(row[10])
				speciess.append(row[7])
			if row[7]=='H&zeta;':
				linetype.append('H_zeta')
				transitions.append(row[10])
				speciess.append(row[7])
		else:
			if fline in fSAC:	
				linetype.append('SAC')
				transitions.append(row[10])
				speciess.append(row[7])
			else:
				if fline in fStrong:
					linetype.append('strong')
					transitions.append(row[10])
					speciess.append(row[7])
				else:
					linetype.append('mol')
					transitions.append(row[10])
					speciess.append(row[7])
		try: fline=float(row[15])
		except:	continue
		if fUnitAF=='GHz': 	flines.append(1000*fline)
		if fUnitAF=='MHz': 	flines.append(fline)
		if fUnitAF=='kHz': 	flines.append(fline/1000.)
		ffits.append(f)
		heights.append(float(row[3]))
		wids.append(float(row[4]))
		if row[16]=='recomb':
			if row[13]=='H&alpha;':
				linetype.append('H_alpha')
			if row[13]=='He&alpha;':
				linetype.append('He_alpha')
			if row[13]=='He&beta;':
				linetype.append('He_beta')
			if row[13]=='H&beta;':
				linetype.append('H_beta')
			if row[13]=='H&gamma;':
				linetype.append('H_gamma')
			if row[13]=='H&delta;':
				linetype.append('H_delta')
			if row[13]=='H&epsilon;':
				linetype.append('H_epsilon')
			if row[13]=='H&zeta;':
				linetype.append('H_zeta')
		else:
			if fline in fSAC:	
				linetype.append('SAC')
				transitions.append(row[16])
				speciess.append(row[13])
			else:
				if fline in fStrong:
					linetype.append('strong')
					transitions.append(row[16])	
					speciess.append(row[13])
				else:
					linetype.append('mol')
					transitions.append(row[16])
					speciess.append(row[13])
	return np.asarray(ffits), np.asarray(flines),np.asarray(linetype),np.asarray(heights), np.asarray(wids),np.asarray(transitions),np.asarray(speciess)


def MASTARslicer(crs,hs,ws,ss,trs,ts,ls,fs,carefulKill='False'):
	crn,hn,wn,sn,trn,tn,ln=[],[],[],[],[],[],[]
	for ii in range(len(crs)):
		if crs[ii] in fs:
			if carefulKill =='False':
				continue
			if carefulKill =='True':
				
				linetype =str(ts[ii])
				if linetype[:2]=='H_' or linetype[:3]=='He_' or linetype =='wing' or linetype =='unidentified':
					continue		
				else:  
					crn.append(crs[ii])
					hn.append(hs[ii])
					wn.append(ws[ii])
					sn.append(ss[ii])
					trn.append(trs[ii])
					tn.append(ts[ii])
					ln.append(ls[ii])
		else:  
			crn.append(crs[ii])
			hn.append(hs[ii])
			wn.append(ws[ii])
			sn.append(ss[ii])
			trn.append(trs[ii])
			tn.append(ts[ii])
			ln.append(ls[ii])
	return np.asarray(crn), np.asarray(hn), np.asarray(wn), np.asarray(sn), np.asarray(trn), np.asarray(tn), np.asarray(ln)

def absmax(ar):
	return np.argmax(np.abs(ar)),ar[np.argmax(np.abs(ar))]

def wingcri1(lf,crV, wV, h,crBig,hBig,wBig,typesBig,vshift):
	crVBig = vshift-(crBig-lf)/lf * c 
	wingmask = abs(crVBig - crV) < wV
	absmask = abs(hBig[wingmask])< 0.22 * abs(h)
	return crVBig[wingmask][absmask],hBig[wingmask][absmask],crBig[wingmask][absmask],wBig[wingmask][absmask],c/lf*wBig[wingmask][absmask],typesBig[wingmask][absmask]

def wingcri2(lf,crV, wV, h,crBig,hBig,wBig,typesBig,vshift):
	crVBig = vshift-(crBig-lf)/lf * c 
	maseFreqs,maserSteerClear=VE.RelevantMaserFreqs('dum')
	if (abs((lf- maseFreqs)/lf*c)).min()<8:
		nowwingmask = abs(crVBig - crV) <maserSteerClear
		absmask = abs(hBig[nowwingmask])< 0.5 * abs(h)
		return crVBig[nowwingmask][absmask],hBig[nowwingmask][absmask],crBig[nowwingmask][absmask],wBig[nowwingmask][absmask],c/lf*wBig[nowwingmask][absmask],typesBig[nowwingmask][absmask]
	else:
		wingmask = abs(crVBig - crV) <1.5* wV
		wingmask1 = abs(crVBig[wingmask] - crV) >1.* wV
		absmask = abs(hBig[wingmask][wingmask1])< 0.05 * abs(h)
		return crVBig[wingmask][wingmask1][absmask],hBig[wingmask][wingmask1][absmask],crBig[wingmask][wingmask1][absmask],wBig[wingmask][wingmask1][absmask],c/lf*wBig[wingmask][wingmask1][absmask],typesBig[wingmask][wingmask1][absmask] 


def compCri(lf,crV, wV, h,crBig,hBig,wBig,typesBig,vshift):
	crVBig = vshift-(crBig-lf)/lf * c 
	wingmask = abs(crVBig - crV) < wV
	absmask = abs(hBig[wingmask])> 0.2* abs(h)
	typemask = typesBig[wingmask][absmask] == 'unidentified'
	return crVBig[wingmask][absmask][typemask],hBig[wingmask][absmask][typemask],crBig[wingmask][absmask][typemask],wBig[wingmask][absmask][typemask],wBig[wingmask][absmask][typemask] * c/lf,typesBig[wingmask][absmask][typemask] 

def wingArrayFill(lf,ls,ltr, wingscrV,wingsh,wingscr,wingsw,wingswV,wingsTypes,lineCR,lspec,ltrans,crVar,har,wVar,crar,ltype,wingOr2Comp='wing'):
	for ii in range(len(wingsh)):
		if wingscr[ii] in crar:
			nnn = np.argmin(abs(crar-wingscr[ii]))
			if wingOr2Comp =='wing':
				ltype[nnn]=wingOr2Comp
		else:
			ltype=np.append(ltype,wingOr2Comp)
			lineCR=np.append(lineCR,lf)
			lspec=np.append(lspec,ls)
			ltrans=np.append(ltrans,ltr)
			crVar=np.append(crVar,wingscrV[ii])
			har=np.append(har,wingsh[ii])
			wVar=np.append(wVar,wingswV[ii])
			crar=np.append(crar,wingscr[ii])
	return np.asarray(lineCR),np.asarray(lspec),np.asarray(ltrans),np.asarray(crVar),np.asarray(har),np.asarray(wVar),np.asarray(crar),np.asarray(ltype)




def wingcheck(place,lineCR,lspec,ltrans,crVar,har,wVar,crar,ltype,crBig,hBig,wBig,typesBig,vshift):
	lenOrig=len(crar)
	dumHar=har
	for dumi in range(lenOrig):
		mm,maxh = absmax(dumHar)
		nn =  np.argmin(abs(har - maxh))
		dumnn = np.argmin(dumHar-maxh)
		dumHar =np.delete(dumHar,mm)
		if ltype[nn]=='wing':
			continue
		dumtype = str(ltype[nn])
		if dumtype[len(dumtype)-6:]=='Recomb':
			dumtype = dumtype[len(dumtype)-6:]
		
		maseFreqs,maserSteerClear=VE.RelevantMaserFreqs('dum')
		if VE.bestcheck(crVar[nn],place,dumtype,har[nn]) =='True' or (abs((lineCR[nn]- maseFreqs)/lineCR[nn]*c)).min()<8:
			lf = lineCR[nn]
			crV = crVar[nn]
			wV = wVar[nn]
			h = har[nn]
			ls=lspec[nn]
			ltr=ltrans[nn]
			

			wingscrV,wingsh,wingscr,wingsw,wingswV,wingsTypes=wingcri1(lf,crV, wV, h,crBig,hBig,wBig,typesBig,vshift)
			if len(wingscrV) > 0:
				lineCR,lspec,ltrans,crVar,har,wVar,crar,ltype = wingArrayFill(lf,ls,ltr, wingscrV,wingsh,wingscr,wingsw,wingswV,wingsTypes,lineCR,lspec,ltrans,crVar,har,wVar,crar,ltype,wingOr2Comp='wing')
			wingscrV,wingsh,wingscr,wingsw,wingswV,wingsTypes=wingcri2(lf,crV, wV, h,crBig,hBig,wBig,typesBig,vshift)
			if len(wingscrV) > 0:
				lineCR,lspec,ltrans,crVar,har,wVar,crar,ltype = wingArrayFill(lf,ls,ltr, wingscrV,wingsh,wingscr,wingsw,wingswV,wingsTypes,lineCR,lspec,ltrans,crVar,har,wVar,crar,ltype,wingOr2Comp='wing')

			comp2scrV,comp2sh,comp2scr,comp2sw,comp2swV,comp2sTypes=compCri(lf,crV, wV, h,crBig,hBig,wBig,typesBig,vshift)
			if len(comp2swV) > 0:
				lineCR,lspec,ltrans,crVar,har,wVar,crar,ltype = wingArrayFill(lf,ls,ltr,comp2scrV,comp2sh,comp2scr,comp2sw,comp2swV,comp2sTypes,lineCR,lspec,ltrans,crVar,har,wVar,crar,ltype,wingOr2Comp='Comp2/Unidentified')
			
			
			
	for ii in range(len(ltype)):
		if ltype[ii] == 'wi':
			ltype[ii] = 'wing'
		if ltype[ii][:3] == 'win':
			ltype[ii] = 'wing'
	return np.asarray(lineCR),np.asarray(lspec),np.asarray(ltrans),np.asarray(crVar),np.asarray(har),np.asarray(wVar),np.asarray(crar),np.asarray(ltype)


def recombFinder(reg, crs,lines,types,hs,wids,transitions,vshift,outwrite='False',plotQ = 'False',o='None'):
	Vs = vshift + c*(lines-crs)/lines
	widsV=wids/crs * c
	def recombKeep(Rcrs,Rvels,Rhs,Rwids,Rlines,Rtrans,rr):
		m=Rhs>0
		med = np.median(Rhs[m])
		m2=Rhs[m]>0.3 * med		
		m3 = Rhs[m][m2]<2. * med
		Rwid2=Rwids[m][m2][m3]
		widmed=np.median(Rwid2)
		m4 = Rwids[m][m2][m3]<2.5*widmed
		m5 = Rwids[m][m2][m3][m4]>0.5*widmed
		Rcrs=Rcrs[m][m2][m3][m4][m5]
		Rvels=Rvels[m][m2][m3][m4][m5]
		Rhs=Rhs[m][m2][m3][m4][m5]
		Rwids=Rwid2[m4][m5]
		Rlines = Rlines[m][m2][m3][m4][m5]
		Rtrans = Rtrans[m][m2][m3][m4][m5]
		dellist=[]
		for dumidumi in range(len(Rtrans)):
			if VE.bestcheck(Rvels[dumidumi],rr,'Recomb',Rhs[dumidumi]) != 'True':
				dellist.append(dumidumi)
		Rcrs=np.delete(Rcrs,dellist)
		Rvels=np.delete(Rvels,dellist)
		Rhs=np.delete(Rhs,dellist)
		Rwids=np.delete(Rwids,dellist)
		Rlines = np.delete(Rlines,dellist)
		Rtrans= np.delete(Rtrans,dellist)
		return Rcrs, Rvels, Rhs, Rwids, Rlines, Rtrans

	Halph=types=='H_alpha'
	Hbet=types=='H_beta'
	Hgam=types=='H_gamma'
	Hdel=types=='H_delta'
	Healph=types=='He_alpha'
	Hebet=types=='He_beta'
	Heps=types=='H_epsilon'
	Hzet=types=='H_zeta'

	HalphCrs, HalphVs, HalphHs,HalphWids,HalphLines,HalphTrans= recombKeep(crs[Halph],Vs[Halph],hs[Halph],widsV[Halph],lines[Halph],transitions[Halph],reg)
	HealphCrs, HealphVs,  HealphHs,HealphWids,HealphLines,HealphTrans = recombKeep(crs[Healph],Vs[Healph],hs[Healph],widsV[Healph],lines[Healph],transitions[Healph],reg)
	HebetCrs, HebetVs,  HebetHs,HebetWids,HebetLines,HebetTrans = recombKeep(crs[Hebet],Vs[Hebet],hs[Hebet],widsV[Hebet],lines[Hebet],transitions[Hebet],reg)
	HbetCrs, HbetVs, HbetHs,HbetWids,HbetLines,HbetTrans = recombKeep(crs[Hbet],Vs[Hbet],hs[Hbet],widsV[Hbet],lines[Hbet],transitions[Hbet],reg)
	HgamCrs,HgamVs,  HgamHs,HgamWids,HgamLines,HgamTrans = recombKeep(crs[Hgam],Vs[Hgam],hs[Hgam],widsV[Hgam],lines[Hgam],transitions[Hgam],reg)
	HdelCrs,HdelVs,  HdelHs,HdelWids,HdelLines,HdelTrans = recombKeep(crs[Hdel],Vs[Hdel],hs[Hdel],widsV[Hdel],lines[Hdel],transitions[Hdel],reg)
	medianDelta = np.median(HdelHs)
	HepsCrs,HepsVs,  HepsHs,HepsWids,HepsLines,HepsTrans = recombKeep(crs[Heps],Vs[Heps],hs[Heps],widsV[Heps],lines[Heps],transitions[Heps],reg)
	medianepsta = np.median(HepsHs)

	Hcrs = np.concatenate((HalphCrs,HbetCrs,HgamCrs,HdelCrs,HepsCrs))
	Hhs = np.concatenate((HalphHs,HbetHs,HgamHs,HdelHs,HepsHs))
	Hwids = np.concatenate((HalphWids,HbetWids,HgamWids,HdelWids,HepsWids))
	HVs = np.concatenate((HalphVs,HbetVs,HgamVs,HdelVs,HepsVs))
	Hlines=np.concatenate((HalphLines,HbetLines,HgamLines,HdelLines,HepsLines))
	Htrans = np.concatenate((HalphTrans,HbetTrans,HgamTrans,HdelTrans,HepsTrans))
	ALLtypes =[]
	ALLspecs=[]
	for ii in range(len(Htrans)):
		ALLtypes.append('H Recomb')
		ALLspecs.append('Hydrogen')
	
	
	ALLcrs = np.concatenate((Hcrs,HealphCrs,HebetCrs))
	ALLhs = np.concatenate((Hhs,HealphHs,HebetHs))
	ALLwids = np.concatenate((Hwids,HealphWids,HebetWids))
	ALLVs = np.concatenate((HVs,HealphVs,HebetVs))
	ALLlines=np.concatenate((Hlines,HealphLines,HebetLines))
	ALLtrans = np.concatenate((Htrans,HealphTrans,HebetTrans))
	for ii in range(len(ALLlines)-len(Hlines)):
		ALLtypes.append('He Recomb')
		ALLspecs.append('Helium')
	ALLtypes = np.asarray(ALLtypes)
	ALLlines,ALLspecs,ALLtrans,ALLVs,ALLhs,ALLwids,ALLcrs,ALLtypes=wingcheck(reg,ALLlines,ALLtypes,ALLtrans,ALLVs,ALLhs,ALLwids,ALLcrs,ALLtypes,crs,hs,wids,types,vshift)

	Hwid_sd= np.sqrt(sum((Hwids-np.mean(Hwids))**2)/(len(Hwids)-1))
	WidMask = abs(Hwids - np.median(Hwids)) < 2 *Hwid_sd
	

	HVs = np.concatenate((HalphVs,HbetVs,HgamVs,HdelVs,HepsVs))
	HVs=HVs[WidMask]
	HV_sd= np.sqrt(sum((HVs-np.mean(HVs))**2)/(len(HVs)-1))
	VMask = abs(HVs - np.median(HVs)) < 2 *HV_sd

	Hwids=Hwids[WidMask][VMask]
	Hwid_sd= np.sqrt(sum((Hwids-np.mean(Hwids))**2)/(len(Hwids)-1))
	HVs=HVs[VMask]
	HV_sd= np.sqrt(sum((HVs-np.mean(HVs))**2)/(len(HVs)-1))
	Hcrs = Hcrs[WidMask][VMask]

	HeGoodV = abs(HealphVs- np.mean(HVs)) < 2 * HV_sd
	HeGoodW = abs(HealphWids[HeGoodV] - np.median(Hwids)) < 2 *Hwid_sd
	HeGoods = HealphCrs[HeGoodV][HeGoodW]
	recomb_good = np.concatenate((Hcrs,HeGoods))
	
	if len(HdelHs) > 3:

		HzetHs=hs[Hzet]
		HzetHeightMask = abs(HzetHs - 0.5 * medianDelta) < (0.2 * medianDelta)
		HzetVs=Vs[Hzet][HzetHeightMask]
		HzetVsMask = abs(HzetVs - np.median(HVs)) < 4 * HV_sd
		HzetWids=widsV[Hzet][HzetHeightMask][HzetVsMask]
		HzetHs=hs[Hzet][HzetHeightMask][HzetVsMask]
		HzetCrs=crs[Hzet][HzetHeightMask][HzetVsMask]
		HzetVs= Vs[Hzet][HzetHeightMask][HzetVsMask]
		HzetLines=lines[Hzet][HzetHeightMask][HzetVsMask]
		HzetTrans=transitions[Hzet][HzetHeightMask][HzetVsMask]
		if len(HzetCrs)>0:
			ALLcrs = np.concatenate((ALLcrs,HzetCrs))
			ALLhs = np.concatenate((ALLhs,HzetHs))
			ALLwids = np.concatenate((ALLwids,HzetWids))
			ALLVs = np.concatenate((ALLVs,HzetVs))
			ALLlines=np.concatenate((ALLlines,HzetLines))
			ALLtrans = np.concatenate((ALLtrans,HzetTrans))
			ALLtypes = list(ALLtypes)
			ALLspecs = list(ALLspecs)
			for ii in range(len(HzetCrs)):
				ALLtypes.append('H Recomb')
				ALLspecs.append('H Recomb')
			ALLtypes = np.asarray(ALLtypes)
			ALLspecs = np.asarray(ALLspecs)
	mastAR = zip(ALLlines,ALLspecs,ALLtrans,ALLVs,ALLhs,ALLwids,ALLcrs,ALLtypes)
	
	if outwrite=='True':
		fileWriter=csv.writer(o,delimiter=':')
		fileWriter.writerow([' ', '', ' ', 'Gaussian Fit','',' ',' '])
		fileWriter.writerow(['Line Rest Frequency (MHz)', 'Species', 'Transition', 'Velocity (km/s)','Height','Width (km/s)','Line Center (MHz)','Line Type'])
		for i in range(len(mastAR)):	
			if mastAR[i][3] in HVs and mastAR[i][5] in Hwids:
				fileWriter.writerow(mastAR[i])
			else:
				if mastAR[i][1]=='He Recomb':	
					fileWriter.writerow(mastAR[i])
				else: 
					dumar =[]
					dumar.extend(mastAR[i])
					dumar.append(['Outlier'])
					fileWriter.writerow(dumar)
		fileWriter.writerow([' ', '', 'H recomb Mean', np.mean(HVs),' ',np.mean(Hwids),' '])
		fileWriter.writerow([' ', '', 'H recomb Median',np.median(HVs),' ',np.median(Hwids),' '])
		fileWriter.writerow([' ', '', 'H recomb Standard Deviation', HV_sd,' ',Hwid_sd,' '])
	return ALLcrs

def SACfinder(reg,fs,hs,ws,mols,trans,types,rfs,vshift,Vrange=[0,0],outwrite='False',plotQ = 'True',r='N',o='None'):
	### This handles line of sight absorption components.  NOT emission.
	try: strongFile,SACfile,fUnitCSV = VE.strongAndSAC(reg)
	except: print "Didn't get anything from VelocityEditFile.py"
	try: fSAC,IDSAC,transSAC = lineGrab(SACfile,fUnitCSV)
	except: 
		print "Could not get data from SACfile.  May not have set a file name.  Proceeding."
		fSAC,IDSAC,transSAC = np.asarray([]),np.asarray([]),np.asarray([])
	
	if Vrange == [0,0]:	
	    if vmain =='dum':
		try:	Vrange = VE.SACvRange(reg)
		except:	
			print "In SACfinder, NEED TO INPUT a valid velocity range for the foreground absorption components.  Otherwise, edit the if statements at start of function to include the region you're targeting"
			return 
	    else:  Vrange,dummyparam = VE.SACvRange(reg)
	print "Began SAC finding!!! Region = ",r
			
	if outwrite=='True':
		fileWriter=csv.writer(o,delimiter=':')
		fileWriter.writerow([' ', '', ' ', 'Gaussian Fit','',' ',' '])
		fileWriter.writerow(['Line Rest Frequency (MHz)', 'Species', 'Transition', 'Velocity (km/s)','Height','Width (km/s)','Line Center (MHz)','Line Type'])
	iiii=0
	SAC_included =[]
	for lf in fSAC:
		vel = vshift-(fs-lf)/lf * c 
		mask1 = vel < Vrange[1]
		mask2 = vel[mask1] > Vrange[0]
		mask3 = hs[mask1][mask2] < 0.0			##Would not see SAC emission if it existed.  
		velComps = vel[mask1][mask2][mask3]
		if len(velComps)<1:
			iiii+=1
			continue
		transs= trans[mask1][mask2][mask3]
		hcomps=hs[mask1][mask2][mask3]
		wcomps= ws[mask1][mask2][mask3]
		fcomps= fs[mask1][mask2][mask3]
		wcompsV=wcomps/lf * c
		lrfs=[]
		transss=[]	
		typess=[]
		spcs=[]
		for jjj in range(len(velComps)):
			lrfs.append(lf)
			typess.append('SAC')
			transss.append(transSAC[iiii])
			spcs.append(IDSAC[iiii])
		lrfs = np.asarray(lrfs)
		transss= np.asarray(transss)
		typess= np.asarray(typess)
		spcs= np.asarray(spcs)
		lrfs,spcs,transss,velComps,hcomps,wcompsV,fcomps,typesss=wingcheck(reg,lrfs,spcs,transss,velComps,hcomps,wcompsV,fcomps,typess,fs,hs,ws,types,vshift)
		OutWriteMastAR=zip(lrfs,spcs,transss,velComps,hcomps,wcompsV,fcomps,typesss)
		velnow =-1000.
		if outwrite=='True':
			for i in range(len(OutWriteMastAR)):
				if velComps[i]!=velnow:
					fileWriter.writerow(OutWriteMastAR[i])
				velnow = velComps[i]
		iiii+=1
		SAC_included.extend(np.unique(fcomps))
	return SAC_included


def strongFinder(reg,fs,hs,ws,mols,trans,types,rfs,vshift,Vrange=[0,0],outwrite='False',plotQ = 'True',r='N',o='None',vmain='dum'):
	print "Began strong line finding!!!", "region = %s " %r
	if Vrange == [0,0]:
	    if vmain =='dum':
		try:	Vrange, vmain = VE.strongVrange(reg)
		except:	
			print "In strongFinder, NEED TO INPUT a valid velocity range for the strong lines.  Otherwise, edit the if statements at start of function to include the region you're targeting"
			return
	    else:  Vrange,dummyparam = VE.strongVrange(reg)
	try: strongFile,SACfile,fUnitCSV = VE.strongAndSAC(reg)
	except: print "Didn't get anything for strongFile & SACfile from VelocityEditFile.py"
	try: fStrong,IDStrong,transStrong = lineGrab(strongFile,fUnitCSV)
	except: 
		print "Could not get data from SACfile.  May not have set a file name.  Proceeding."
		fStrong,IDStrong,transStrong = np.asarray([]),np.asarray([]),np.asarray([])
	if outwrite=='True':
		fileWriter=csv.writer(o,delimiter=':')
		fileWriter.writerow([' ', '', ' ', 'Gaussian Fit','',' ',' '])
		fileWriter.writerow(['Line Rest Frequency (MHz)', 'Species', 'Transition', 'Velocity (km/s)','Height','Width (km/s)','Line Center (MHz)','Line Type'])
	iiii=0
	fs_LSR = fs - vshift/c * fs
	strongs_included=[]
	for lf in fStrong:
		vel = (lf - fs_LSR)/lf * c 
		mask1 = vel < Vrange[1]
		mask2 = vel[mask1] > Vrange[0]
		velComps = vel[mask1][mask2]
		if len(velComps)<1:
			iiii+=1
			continue
		transs= trans[mask1][mask2]
		hcomps=hs[mask1][mask2]
		wcomps= ws[mask1][mask2]
		fcomps= fs[mask1][mask2]
		wcompsV=wcomps/lf * c
		rfcomps=rfs[mask1][mask2]

		mdum= transs!='un'
		theTranss=transs[mdum]
		velTranss=vel[mdum]
		if len(theTranss)>0:
			min64 = np.argmin(abs(velTranss-vmain))
			if rfcomps[min64]!= lf:
				min64 = np.argmin(abs(rfcomps[mdum] - lf))
			theTrans = theTranss[min64]
		else:  theTrans ='dummy'	
		lrfs=[]
		transss=[]	
		typess=[]
		spcs=[]
		for jjj in range(len(velComps)):
			lrfs.append(lf)
			typess.append('strong')
			transss.append(theTrans)
			spcs.append(IDStrong[iiii])
		#print 'transss',transss
	#	print "IM IN strongFinder!! lrfs = ",lrfs
	#	print "		  velcomps = ",velComps
		lrfs = np.asarray(lrfs)
		transss= np.asarray(transss)
		typess= np.asarray(typess)
		spcs= np.asarray(spcs)
		lrfs,spcs,transss,velComps,hcomps,wcompsV,fcomps,typess=wingcheck(reg,lrfs,spcs,transss,velComps,hcomps,wcompsV,fcomps,typess,fs,hs,ws,types,vshift)
		strongs_included.extend(np.unique(fcomps))
		OutWriteMastAR=zip(lrfs,spcs,transss,velComps,hcomps,wcompsV,fcomps,typess)
		velnow =-1000.
	#	print 'len(OutWriteMastAR)', len(OutWriteMastAR)
	#	print 'OutWriteMastAR', OutWriteMastAR
		if outwrite=='True':
			for i in range(len(OutWriteMastAR)):
				if velComps[i]!=velnow:
					#print 'OutWriteMastAR[i]',OutWriteMastAR[i]
					fileWriter.writerow(OutWriteMastAR[i])
				velnow = velComps[i]
		iiii+=1
	return strongs_included

def molID(reg,fs,hs,ws,mols,trans,types,rfs,vshift,outwrite='False',plotQ = 'True',r='N',o='None'):
	## reg = region you're looking at.  fs = frequencies, hs=heights, ws=widths, mols=molecules, trans=transitions.  
	###If your frequency axis is shifted to some rest velocity, vshift = restvelocity
	print "Begin 1-component Line ID-ing!!!", "region = %s " %r
	
	if outwrite=='True':
		fileWriter=csv.writer(o,delimiter=':')
		fileWriter.writerow([' ', '', ' ', 'Gaussian Fit','',' ',' '])
		fileWriter.writerow(['Line Rest Frequency (MHz)', 'Species', 'Transition', 'Velocity (km/s)','Height','Width (km/s)','Line Center (MHz)','Line Type'])
	iiii=0
	mask = types == 'mol'
	fsOrig=fs
	hsOrig=hs
	wsOrig=ws
#	print "IN MOL.  RFS to start with!:",rfs
#	print "fs to start with:",fs
#	print "types to start with",types
	typesOrig=types
	fs = fs[mask]
	hs = hs[mask]
	ws = ws[mask]
	mols = mols[mask]
	trans = trans[mask]
	types = types[mask]
	rfs = rfs[mask]
	fs_LSR = fs - vshift/c * fs
	fskill = []
#	print 'rfs in mol', rfs
#	print "IN MOL.  WANT TO KNOW WHAT'S IN rfs", rfs
	for lf in rfs:
		if lf ==0:
			continue
		#print 'iiii',iiii
		vel = (lf - fs_LSR[iiii])/lf * c 
		lrfs,spcs,transss,velComps,hcomps,wcompsV,fcomps,typess=wingcheck(reg,np.asarray([lf]),np.asarray([mols[iiii]]),np.asarray([trans[iiii]]),np.asarray([vel]),np.asarray([hs[iiii]]),np.asarray([ws[iiii]])/lf * c,np.asarray([fs[iiii]]),np.asarray([types[iiii]]),fsOrig,hsOrig,wsOrig,typesOrig,vshift)
		OutWriteMastAR=zip(lrfs,spcs,transss,velComps,hcomps,wcompsV,fcomps,typess)
		velnow =-1000.
		if outwrite=='True':
			for i in range(len(OutWriteMastAR)):
				if velComps[i]!=velnow and fcomps[i] not in fskill:
					#print 'OutWriteMastAR[i]',OutWriteMastAR[i]
					fileWriter.writerow(OutWriteMastAR[i])
				velnow = velComps[i]
		fskill.extend(np.unique(fcomps))
		iiii+=1
	return fskill




def otherU(reg,fs,hs,ws,mols,trans,types,rfs,vshift,outwrite='False',plotQ = 'True',r='N',o = 'None'):
	print "Begin the remaining lines!!!", "region = %s " %r
	if outwrite=='True':
		fileWriter=csv.writer(o,delimiter=':')
		fileWriter.writerow([' ', '', ' ', 'Gaussian Fit','',' ',' '])
		fileWriter.writerow(['Line Rest Frequency (MHz)', 'Species', 'Transition', 'Velocity (km/s)','Height','Width (km/s)','Line Center (MHz)','Line Type'])
	fs_LSR = fs - vshift/c * fs
	iiii=0
	for lf in rfs:
		#print 'lf!',lf, 'Species',mols[iiii]
		#print 'iiii',iiii
		if outwrite == 'True':
			if lf!=0:
			 	vel = (lf - fs_LSR[iiii])/lf * c 
				fileWriter.writerow([lf, mols[iiii],trans[iiii],vel,hs[iiii],ws[iiii]/fs[iiii] * c,fs[iiii],'unidentified'])
			else:
				fileWriter.writerow([lf, mols[iiii],trans[iiii],'---',hs[iiii],ws[iiii]/fs[iiii] * c,fs[iiii],'unidentified'])
		iiii+=1
	return fs



###############################################################################################
#######################    NOW WE RUN THROUGH OUR DATA    #####################################
REGS = VE.REGIONS('i')

for r in REGS:
	print '################################################################################'
	print 'region = %s !!!~!!!' % r
	
	allFitsFile,ALLFITSUnits= VE.ALLFITSin(r)
	VShift=VE.velshift(r)
	CRs,Lines,Types,Hs,Wids,TRans,Specs=ArrayFiller(r,allFitsFile,fUnitAF=ALLFITSUnits)
	print "len of ArrayFiller arrays",len(CRs), len(Wids), len(Types)
	
	outfile=open('velocity_ALL_%s.csv' % r,'wt')
	print "\n GOING INTO RECOMBS"
	RecombKill= recombFinder(r,CRs,Lines,Types,Hs,Wids,TRans,VShift,outwrite='True',plotQ = 'False',o=outfile)
	############ NOTE:  THIS FINDS THE USUAL LINE PARAMETERS AND USES THOSE TO DETERMINE IF SOME LINES ARE NOT CONSISTENT WITH RECOMBINATION LINES.  IT ONLY CONSIDERS RECOMBINATION LINE EMISSION (NO ABSORPTION) AND IT WORKS BEST WHEN YOU HAVE A LARGE SAMPLE OF LINES.
	CRs,Hs,Wids,Specs,TRans,Types,Lines = MASTARslicer(CRs,Hs,Wids,Specs,TRans,Types,Lines, RecombKill,carefulKill ='False')

	print "\n HANDLING STRONGLINES"
	vR,vm = VE.strongVrange(r)
	strong_accounted4 = strongFinder(r,CRs,Hs,Wids,Specs,TRans,Types,Lines,VShift,Vrange=vR,outwrite='True',plotQ = 'True',r='%s' % r,o=outfile,vmain=vm)
	############ NOTE:  THIS FINDS LINES THAT IVE ISOLATED AS "MULTIPLE VELOCITY LINES", WHICH ARE STRONG AND IN MULTIPLE VELOCITY COMPONENTS.  CAN INPUT vshift, vrange, and vmain to change the velocity it's shifted by and the range of velocities in which it will look line radiation, and the .
	CRs,Hs,Wids,Specs,TRans,Types,Lines = MASTARslicer(CRs,Hs,Wids,Specs,TRans,Types,Lines, strong_accounted4)

	print "\n HANDLING FOREGROUND ABSORPTION"
	vR = VE.SACvRange(r)
	SAC_accounted4 = SACfinder(r,CRs,Hs,Wids,Specs,TRans,Types,Lines,VShift,Vrange=vR,outwrite='True',plotQ = 'True',r='%s' % r,o=outfile)
	############ NOTE:  THIS FINDS LINES THAT IVE ISOLATED AS BEING PRESENT IN THE SPIRAL ARM CLOUDS.  CAN INPUT vshift and vrange to change the velocity it's shifted by and the range of velocities in which it will look for line-of-sight absorption components.
	CRs,Hs,Wids,Specs,TRans,Types,Lines = MASTARslicer(CRs,Hs,Wids,Specs,TRans,Types,Lines, SAC_accounted4)

	print "\n GOING INTO MOLS"
	molecsIDd = molID(r,CRs,Hs,Wids,Specs,TRans,Types,Lines,VShift,outwrite='True',plotQ = 'True',r='%s' % r,o=outfile)
	############ NOTE:  THIS FINDS THE VELOCITY OF THE REMAINING MOLECULAR LINES.  DOES THIS LAST BECAUSE WE MAY BE LESS CONFIDENT IN THESE LINE ID's.  
	CRs,Hs,Wids,Specs,TRans,Types,Lines = MASTARslicer(CRs,Hs,Wids,Specs,TRans,Types,Lines, molecsIDd)

	print "\n GOING INTO UNIDENTIFIEDS"
	UIDd = otherU(r,CRs,Hs,Wids,Specs,TRans,Types,Lines,VShift,outwrite='True',plotQ = 'True',r='%s' % r,o=outfile)
	############ NOTE:  THIS PUTS THE REMAINING LINES INTO OUR OUTPUT CSV FILE.  IT WILL INCLUDE THE INFORMATION ABOUT THE TENTATIVE/INCONCLUSIVE TRANSITIONS THAT ARE NEARBY.  














