import numpy as np, pylab as pl, csv
 
c=2.99792E5
regions =('N','h','K4','L','LMH','blank','M')
MethMaserFreq = 36169.29		#in MHz
#coolbins=[0,1,4]
coolbins=[0,1,2,3,4]

##NOTE:  I WILL NEED TO REFIT ANY LINES WITH BINCODE = 5 IN ALLFITS BEFORE RUNNING THIS!!!!

#regions =('K4','K6','L')

regions =['K6','K4','L']

#####################   DEFINE A BUNCH OF FUNCTIONS.  SEE COMMENTS IN the functions Bestcheck, 
def freqGrab(infile):
	tt=csv.reader(open(infile,'rt'),delimiter=':')
	fs=[]
	for row in tt:
		try: 	f=float(row[0])
		except:	continue
		fs.append(f)		
	return np.asarray(fs)
	
def IDGrab(infile):
	tt=csv.reader(open(infile,'rt'),delimiter=':')
	IDs=[]
	for row in tt:
		try: 	f=float(row[0])
		except:	continue
		IDs.append(row[1])		
	return np.asarray(IDs)
def transGrab(infile):
	tt=csv.reader(open(infile,'rt'),delimiter=':')
	IDs=[]
	for row in tt:
		try: 	f=float(row[0])
		except:	continue
		IDs.append(row[2])		
	return np.asarray(IDs)
				
def ArrayFiller(editfile):
	t=csv.reader(open(editfile,'rt'),delimiter=':')
	fvel = freqGrab('multipleVel.csv')
	IDvel=IDGrab('multipleVel.csv')
	fSAC = freqGrab('SAClines.csv')
	Umask= IDvel=='U'
	fvelU=fvel[Umask]
	ffits=[]
	flines=[]
	linetype=[]
	heights=[]
	wids=[]
	transitions=[]
	speciess=[]
	iii = 0
	for row in t:	
		try: f=float(row[1])	#I had it in GHz.  Currently plotting in MHz.
		except: continue
		if int(row[6]) not in coolbins:
			continue
		ffits.append(f)
		heights.append(float(row[3]))
		wids.append(float(row[4]))
		try: fline=float(row[9])	#I had it in GHz.  Currently plotting in MHz.
		except: 
			iii+=1
			if np.min(abs(fvelU-f))<.002:
				linetype.append('velU')	
				flines.append(1000*fvelU[np.argmin(abs(fvelU-f))])
				transitions.append(row[10])
				speciess.append(row[7])
			else:
				linetype.append('unidentified')
				flines.append(0.0)
				transitions.append('un')
				speciess.append('un')
			continue
		flines.append(1000*fline)
		if row[11]=='Recomb':
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
				if fline in fvel:
					linetype.append('vel')
					transitions.append(row[10])
					speciess.append(row[7])
				else:
					linetype.append('mol')
					transitions.append(row[10])
					speciess.append(row[7])
		try: fline=float(row[14])
		except:	continue
		flines.append(1000*fline)
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
				transitions.append(row[15])
				speciess.append(row[12])
			else:
				if fline in fvel:
					linetype.append('vel')
					transitions.append(row[15])	
					speciess.append(row[12])
				else:
					linetype.append('mol')
					transitions.append(row[15])
					speciess.append(row[12])
	return np.asarray(ffits), np.asarray(flines),np.asarray(linetype),np.asarray(heights), np.asarray(wids),np.asarray(transitions),np.asarray(speciess)

def MASTARslicer(crs,hs,ws,ss,trs,ts,ls,fs,carefulKill='False'):
	crn=[]
	hn=[]
	wn=[]
	sn=[]
	trn=[]
	tn=[]
	ln=[]
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

##########  Either edit bestcheck to include your region and line types or input a velcr and a velwid  ###########
def bestcheck(vel,region='K6',linetype='dum',height=0,velcr='dum',velwid='dum'):
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
		if linetype =='mol':		## The Mol's seem to be all over the place... what's the issue?  
			return 'True'
		if linetype == 'mV': 
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
		if linetype == 'mV': 
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
		if linetype == 'mV': 
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
		if linetype == 'mV': 
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
		if linetype == 'mV': 
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
		if linetype == 'mV': 
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

def absmax(ar):
	return np.argmax(np.abs(ar)),ar[np.argmax(np.abs(ar))]

def wingcri1(lf,crV, wV, h,crBig,hBig,wBig,typesBig,vshift=64):
	crVBig = vshift-(crBig-lf)/lf * c 
	wingmask = abs(crVBig - crV) < wV
	absmask = abs(hBig[wingmask])< 0.22 * abs(h)
	return crVBig[wingmask][absmask],hBig[wingmask][absmask],crBig[wingmask][absmask],wBig[wingmask][absmask],c/lf*wBig[wingmask][absmask],typesBig[wingmask][absmask]

def wingcri2(lf,crV, wV, h,crBig,hBig,wBig,typesBig):
	crVBig = 64.-(crBig-lf)/lf * c 
	if abs(lf- MethMaserFreq)<0.5:
		nowwingmask = abs(crVBig - crV) <50.
		absmask = abs(hBig[nowwingmask])< 0.5 * abs(h)
		return crVBig[nowwingmask][absmask],hBig[nowwingmask][absmask],crBig[nowwingmask][absmask],wBig[nowwingmask][absmask],c/lf*wBig[nowwingmask][absmask],typesBig[nowwingmask][absmask]
	else:
		wingmask = abs(crVBig - crV) <1.5* wV
		wingmask1 = abs(crVBig[wingmask] - crV) >1.* wV
		absmask = abs(hBig[wingmask][wingmask1])< 0.05 * abs(h)
		return crVBig[wingmask][wingmask1][absmask],hBig[wingmask][wingmask1][absmask],crBig[wingmask][wingmask1][absmask],wBig[wingmask][wingmask1][absmask],c/lf*wBig[wingmask][wingmask1][absmask],typesBig[wingmask][wingmask1][absmask] 


def compCri(lf,crV, wV, h,crBig,hBig,wBig,typesBig,vshift=64):
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
	return lineCR,lspec,ltrans,crVar,har,wVar,crar,ltype

def wingcheck(place,lineCR,lspec,ltrans,crVar,har,wVar,crar,ltype,crBig,hBig,wBig,typesBig):
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
		if bestcheck(crVar[nn],place,dumtype,har[nn]) =='True' or abs(lineCR[nn] - MethMaserFreq)<1.:
	#		print "  "
	#		print "  "
	#		print "  "
	#		print "Bestcheck = true for mm, maxh",mm,maxh
			lf = lineCR[nn]
			crV = crVar[nn]
			wV = wVar[nn]
			h = har[nn]
			ls=lspec[nn]
			ltr=ltrans[nn]
			
	#		print "crV =",crV

			wingscrV,wingsh,wingscr,wingsw,wingswV,wingsTypes=wingcri1(lf,crV, wV, h,crBig,hBig,wBig,typesBig)
			if len(wingscrV) > 0:
	#			print "got things from wingcriterion 1"
	#			print "wingscrV",wingscrV
	#			print "wingsh",wingsh
	#			print 'len(wingscrV) > 0!',  len(wingscrV)
	#			print "Before wingArrayFill, len(lineCR) =",len(lineCR)
				lineCR,lspec,ltrans,crVar,har,wVar,crar,ltype = wingArrayFill(lf,ls,ltr, wingscrV,wingsh,wingscr,wingsw,wingswV,wingsTypes,lineCR,lspec,ltrans,crVar,har,wVar,crar,ltype,wingOr2Comp='wing')
	#			print 'after criterion 1, len(lineCR) = ',  len(lineCR)
	#			print "lineCR=",lineCR
	#			print "crVar",crVar
			wingscrV,wingsh,wingscr,wingsw,wingswV,wingsTypes=wingcri2(lf,crV, wV, h,crBig,hBig,wBig,typesBig)
			if len(wingscrV) > 0:
	#			print "got things from wingcriterion 2"
	#			print 'len(wingscrV) > 0!',  len(wingscrV)
				lineCR,lspec,ltrans,crVar,har,wVar,crar,ltype = wingArrayFill(lf,ls,ltr, wingscrV,wingsh,wingscr,wingsw,wingswV,wingsTypes,lineCR,lspec,ltrans,crVar,har,wVar,crar,ltype,wingOr2Comp='wing')

			comp2scrV,comp2sh,comp2scr,comp2sw,comp2swV,comp2sTypes=compCri(lf,crV, wV, h,crBig,hBig,wBig,typesBig)
			if len(comp2swV) > 0:
				#print 'len(comp2swV) > 0!',  len(comp2swV)
				lineCR,lspec,ltrans,crVar,har,wVar,crar,ltype = wingArrayFill(lf,ls,ltr,comp2scrV,comp2sh,comp2scr,comp2sw,comp2swV,comp2sTypes,lineCR,lspec,ltrans,crVar,har,wVar,crar,ltype,wingOr2Comp='Comp2/Unidentified')
			
			
			
#	print "ltype before!", ltype
	for ii in range(len(ltype)):
		if ltype[ii] == 'wi':
			ltype[ii] = 'wing'
		if ltype[ii][:3] == 'win':
			ltype[ii] = 'wing'
#	print "ltype now!", ltype
	return np.asarray(lineCR),np.asarray(lspec),np.asarray(ltrans),np.asarray(crVar),np.asarray(har),np.asarray(wVar),np.asarray(crar),np.asarray(ltype)

def recombPlotter(reg, crs,lines,types,hs,wids,transitions,outwrite='False',plotQ = 'True',o='None',vshift=64):
	Vs = vshift + c*(lines-crs)/lines
	widsV=wids/crs * c
	def recombKeep(Rcrs,Rvels,Rhs,Rwids,Rlines,Rtrans,rr):
		m=Rhs>0
		med = np.median(Rhs[m])
		m2=Rhs[m]>0.3 * med		##Dec 24.  I changed this from 0.5 to 0.3 because its cutting out the low resolution recomb lines for L and maybe K6.  
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
			####### Note, you could do "if bestcheck(velcr,velwid)!='True':"
			if bestcheck(Rvels[dumidumi],rr,'Recomb',Rhs[dumidumi]) != 'True':
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
	ALLlines,ALLspecs,ALLtrans,ALLVs,ALLhs,ALLwids,ALLcrs,ALLtypes=wingcheck(reg,ALLlines,ALLtypes,ALLtrans,ALLVs,ALLhs,ALLwids,ALLcrs,ALLtypes,crs,hs,wids,types)

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
	#print "ALLcrs after!!", len(ALLcrs)
	#print "ALLcrs after!!", ALLcrs
	#print 'len(FinalKeep)', len(Hwids)
	#print ' 	'

	#print 'Hwids',Hwids
	#print 'np.mean(Hwids)',np.mean(Hwids)
	#print 'np.mean(Hwids)',np.mean(Hwids)
	#print 'Hwid_sd',Hwid_sd
	#print '	'
	#print 'HVs',HVs
	#print 'np.median(HVs)',np.median(HVs)
	#print 'np.mean(HVs)',np.mean(HVs)
	#print 'HV_sd',HV_sd
	
	if outwrite=='True':
		fileWriter=csv.writer(o,delimiter=':')
		fileWriter.writerow([' ', '', ' ', 'Gaussian FIt','',' ',' '])
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
	#return recomb_good, ALLcrs
	return ALLcrs


def SACfinder(reg,fs,hs,ws,mols,trans,types,rfs,outwrite='False',plotQ = 'True',r='N',o='None',vshift=64,Vrange=[0,0],vmain='dum'):
	### This handles line of sight absorption components.  NOT emission.
	print "Began SAC finding!!! Region = ",r
	SACrfs = 1000*freqGrab('SAClines.csv')
	SACspecs=IDGrab('SAClines.csv')
	SACtrans=transGrab('SAClines.csv')
	if r =='L':
		Vrange=[-120,50]
	if r =='LMH':
		Vrange=[-120,55]
	if r == 'K6' or r == 'N':
		Vrange =[-120,58]
	if r =='K4':
		Vrange=[-20,60]
	if r =='M':
		Vrange=[-120,40]
	if Vrange == [0,0]:
		print "NEED TO INPUT A valid range in velocity for the multiple-velocity lines.  Otherwise, edit the if statements at start of function to include the region you're targeting"
		return 
			
	if outwrite=='True':
		fileWriter=csv.writer(o,delimiter=':')
		fileWriter.writerow([' ', '', ' ', 'Gaussian FIt','',' ',' '])
		fileWriter.writerow(['Line Rest Frequency (MHz)', 'Species', 'Transition', 'Velocity (km/s)','Height','Width (km/s)','Line Center (MHz)','Line Type'])
	iiii=0
	SAC_included =[]
	for lf in SACrfs:
		#print 'lf!',lf, 'Species',SACspecs[iiii]
		#print 'iiii',iiii
		vel = vshift-(fs-lf)/lf * c 
		mask1 = vel < Vrange[1]
		mask2 = vel[mask1] > Vrange[0]
		mask3 = hs[mask1][mask2] < 0.0			##Would not see SAC emission if it existed.  
		#print 'vel[mask1][mask2]',vel[mask1][mask2]
		#print 'fs[mask1][mask2]',fs[mask1][mask2]
		#print 'hs[mask1][mask2]',hs[mask1][mask2]
		velComps = vel[mask1][mask2][mask3]
		#print 'vel[mask1][mask2][mask3]',vel[mask1][mask2][mask3]
		#print 'fs[mask1][mask2][mask3]',fs[mask1][mask2][mask3]
		if len(velComps)<1:
			iiii+=1
			continue
		#print 'velComps!!!', np.unique(velComps)
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
			transss.append(SACtrans[iiii])
			spcs.append(SACspecs[iiii])
		lrfs = np.asarray(lrfs)
		transss= np.asarray(transss)
		typess= np.asarray(typess)
		spcs= np.asarray(spcs)
	#	print "fcomps!!!! =",fcomps
	#	print 'lf', lf, 'species', spcs
	#	print 'Im in mV: len(lrfs) before:', len(lrfs)
		lrfs,spcs,transss,velComps,hcomps,wcompsV,fcomps,typesss=wingcheck(reg,lrfs,spcs,transss,velComps,hcomps,wcompsV,fcomps,typess,fs,hs,ws,types)
	#	print "typessssssssss!!!!! = ",typesss
	#	print 'Im in SAC: len(lrfs) after:', len(lrfs)
		OutWriteMastAR=zip(lrfs,spcs,transss,velComps,hcomps,wcompsV,fcomps,typesss)
		velnow =-1000.
		if outwrite=='True':
			for i in range(len(OutWriteMastAR)):
				if velComps[i]!=velnow:
				#	print 'OutWriteMastAR[i]',OutWriteMastAR[i]
					fileWriter.writerow(OutWriteMastAR[i])
				velnow = velComps[i]
		iiii+=1
		SAC_included.extend(np.unique(fcomps))
	return SAC_included

def velFinder(reg,fs,hs,ws,mols,trans,types,rfs,outwrite='False',plotQ = 'True',r='N',o='None',vshift=64,Vrange=[0,0],vmain='dum'):
	if vmain == 'dum':
		vmain = vshift    ### VMAIN IS THE STRONGEST OR PRIMARY VELOCITY COMPONENT.  
	print "Began 2-component finding!!!", "region = %s " %r
	mVrfs = 1000*freqGrab('multipleVel.csv')
	mVspecs=IDGrab('multipleVel.csv')
	fs_LSR = fs - vshift/c * fs
	if r =='L':
		Vrange=[50,80]
		vmain = 55
	if r =='LMH':
		Vrange=[55,85]
	if r == 'K6' or r == 'N': 
		Vrange =[58,85]
	if r =='K4':
		Vrange =[58,88]
	if r == 'h':
		Vrange =[60,88]
	if r == 'M':
		Vrange =[40,70]
		vmain = 50
	if Vrange == [0,0]:
		print "NEED TO INPUT A valid range in velocity for the multiple-velocity lines.  Otherwise, edit the if statements at start of function to include the region you're targeting"
		return 
			
	if outwrite=='True':
		fileWriter=csv.writer(o,delimiter=':')
		fileWriter.writerow([' ', '', ' ', 'Gaussian FIt','',' ',' '])
		fileWriter.writerow(['Line Rest Frequency (MHz)', 'Species', 'Transition', 'Velocity (km/s)','Height (MJy/sr)','Width (km/s)','Line Center (MHz)','Line Type'])
	iiii=0
	mV_included=[]
	print 'VRANGE!!!!', Vrange
	for lf in mVrfs:
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
			typess.append('mV')
			transss.append(theTrans)
			spcs.append(mVspecs[iiii])
		#print 'transss',transss
	#	print "IM IN MV!! lrfs = ",lrfs
	#	print "		  velcomps = ",velComps
		lrfs = np.asarray(lrfs)
		transss= np.asarray(transss)
		typess= np.asarray(typess)
		spcs= np.asarray(spcs)
		#print 'lf', lf, 'species', spcs
		#print 'Im in mV: len(lrfs) before:', len(lrfs)
		lrfs,spcs,transss,velComps,hcomps,wcompsV,fcomps,typess=wingcheck(reg,lrfs,spcs,transss,velComps,hcomps,wcompsV,fcomps,typess,fs,hs,ws,types)
		mV_included.extend(np.unique(fcomps))
		#print 'Im in mV: len(lrfs) after:', len(lrfs)
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
	return mV_included


def molID(reg,fs,hs,ws,mols,trans,types,rfs,outwrite='False',plotQ = 'True',r='N',o='None',vshift = 64):
	## reg = region you're looking at.  fs = frequencies, hs=heights, ws=widths, mols=molecules, trans=transitions.  
	###If your frequency axis is shifted to some rest velocity, vshift = restvelocity
	print "Begin 1-component Line ID-ing!!!", "region = %s " %r
	
	if outwrite=='True':
		fileWriter=csv.writer(o,delimiter=':')
		fileWriter.writerow([' ', '', ' ', 'Gaussian FIt','',' ',' '])
		fileWriter.writerow(['Line Rest Frequency (MHz)', 'Species', 'Transition', 'Velocity (km/s)','Height (MJy/sr)','Width (km/s)','Line Center (MHz)','Line Type'])
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
		lrfs,spcs,transss,velComps,hcomps,wcompsV,fcomps,typess=wingcheck(reg,np.asarray([lf]),np.asarray([mols[iiii]]),np.asarray([trans[iiii]]),np.asarray([vel]),np.asarray([hs[iiii]]),np.asarray([ws[iiii]])/lf * c,np.asarray([fs[iiii]]),np.asarray([types[iiii]]),fsOrig,hsOrig,wsOrig,typesOrig)
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

def otherU(reg,fs,hs,ws,mols,trans,types,rfs,outwrite='False',plotQ = 'True',r='N',o = 'None',vshift = 64):
	print "Begin the remaining lines!!!", "region = %s " %r
	if outwrite=='True':
		fileWriter=csv.writer(o,delimiter=':')
		fileWriter.writerow([' ', '', ' ', 'Gaussian FIt','',' ',' '])
		fileWriter.writerow(['Line Rest Frequency (MHz)', 'Species', 'Transition', 'Velocity (km/s)','Height (MJy/sr)','Width (km/s)','Line Center (MHz)','Line Type'])
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
for r in regions:
	print '	'
	print 'region = %s !!!~!!!' % r
	CRs,Lines,Types,Hs,Wids,TRans,Specs=ArrayFiller('../ALLFITS_%s_Sept14_4.csv' % r)
	print "len of ArrayFiller arrays",len(CRs), len(Wids), len(Types)
	
	outfile=open('velocity_ALL_Sept14_4_%s.csv' % r,'wt')
	print "\n GOING INTO RECOMBS"
	RecombKill= recombPlotter(r,CRs,Lines,Types,Hs,Wids,TRans,outwrite='True',plotQ = 'False',o=outfile)
	############ NOTE:  THIS FINDS THE USUAL LINE PARAMETERS AND USES THOSE TO DETERMINE IF SOME LINES ARE NOT CONSISTENT WITH RECOMBINATION LINES.  IT ONLY CONSIDERS RECOM LINE EMISSION (NO ABSORPTION) AND IT WORKS BEST WHEN YOU HAVE A LARGE SAMPLE OF LINES.
	CRs,Hs,Wids,Specs,TRans,Types,Lines = MASTARslicer(CRs,Hs,Wids,Specs,TRans,Types,Lines, RecombKill,carefulKill ='False')

	print "\n GOING INTO mV"
	mV_accounted4 = velFinder(r,CRs,Hs,Wids,Specs,TRans,Types,Lines,outwrite='True',plotQ = 'True',r='%s' % r,o=outfile)
	############ NOTE:  THIS FINDS LINES THAT IVE ISOLATED AS "MULTIPLE VELOCITY LINES", WHICH ARE STRONG AND IN MULTIPLE VELOCITY COMPONENTS.  CAN INPUT vshift, vrange, and vmain to change the velocity it's shifted by and the range of velocities in which it will look line radiation, and the .
	CRs,Hs,Wids,Specs,TRans,Types,Lines = MASTARslicer(CRs,Hs,Wids,Specs,TRans,Types,Lines, mV_accounted4)

	print "\n GOING INTO SAC"
	SAC_accounted4 = SACfinder(r,CRs,Hs,Wids,Specs,TRans,Types,Lines,outwrite='True',plotQ = 'True',r='%s' % r,o=outfile)
	############ NOTE:  THIS FINDS LINES THAT IVE ISOLATED AS BEING PRESENT IN THE SPIRAL ARM CLOUDS.  CAN INPUT vshift and vrange to change the velocity it's shifted by and the range of velocities in which it will look for line-of-sight absorption components.
	CRs,Hs,Wids,Specs,TRans,Types,Lines = MASTARslicer(CRs,Hs,Wids,Specs,TRans,Types,Lines, SAC_accounted4)

	print "\n GOING INTO MOLS"
	molecsIDd = molID(r,CRs,Hs,Wids,Specs,TRans,Types,Lines,outwrite='True',plotQ = 'True',r='%s' % r,o=outfile)
	############ NOTE:  THIS FINDS THE VELOCITY OF THE REMAINING MOLECULAR LINES.  DOES THIS LAST BECAUSE WE MAY BE LESS CONFIDENT IN THESE LINE ID's.  
	CRs,Hs,Wids,Specs,TRans,Types,Lines = MASTARslicer(CRs,Hs,Wids,Specs,TRans,Types,Lines, molecsIDd)

	print "\n GOING INTO UNIDENTIFIEDS"
	UIDd = otherU(r,CRs,Hs,Wids,Specs,TRans,Types,Lines,outwrite='True',plotQ = 'True',r='%s' % r,o=outfile)
	############ NOTE:  THIS PUTS THE REMAINING LINES INTO OUR OUTPUT CSV FILE.  IT WILL INCLUDE THE INFORMATION ABOUT THE TENTATIVE/INCONCLUSIVE TRANSITIONS THAT ARE NEARBY.  
