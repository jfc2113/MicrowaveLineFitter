import numpy as np, pylab as pl, csv
from math import *

############### I need to put this in terms of crvel & velrange.  

def lineCompare(infile,outfile,bigAr,crAr,plotQ=False,writeout=True):
	print "aggh!"
	if infile=='None':
	
	    if writeout ==True or writeout == 'True':
		o=open(outfile,'wt')
		fileWriter=csv.writer(o,delimiter=':')
		for m in range(len(bigAr)):
			params=bigAr[m]	
			fileWriter.writerow(params)
	else:
		t=csv.reader(open(infile,'rb'),delimiter=':')
		if writeout ==True or writeout == 'True':
			o=open(outfile,'wt')
			fileWriter=csv.writer(o,delimiter=':')
		flines=[]
		trans=[]
		for row in t:	
			try: f=1000*float(row[2])	#I had it in GHz.  Currently plotting in MHz.
			except: continue
			trans=row[3]
		
			if plotQ==True or plotQ == 'True':
				pl.plot([f,f],[-1*500,500],':',color='r')
			if np.abs(crAr-f).min() < 2:
				IDmask=abs(crAr-f)<2.
				args = np.arange(len(crAr))[IDmask]
				print "but here?"
				for m in args:
					dumar=[]
					dumar.append(f)
					dumar.append(trans)
					print "goodness!"
					params=bigAr[m]		
					##Should be the cr, croffset, height, fwhm, and list membership 
					print "params!", params
					dumar.extend(params)
					print "dumar!",dumar 
					if writeout == True or writeout == 'True':
						fileWriter.writerow(dumar)
			else:
				   if writeout == True or writeout == 'True':
					dumar=[]
					dumar.append(f)
					dumar.append(trans)
					fileWriter.writerow(dumar)
	if writeout== True or writeout == 'True':			
		o.close()


def lineCompareAll(infile,outfile,bigAr,plotQ=False,writeout=True,infile2 = 'None',velcr = 72,delVel=17):
	c = 2.9979E5 #  in km/s	
	## Output will indicate lines that are within delVel km/s of velcr.  
	if infile=='None':
	    if writeout == True or writeout == 'True':
		o=open(outfile,'wt')
		fileWriter=csv.writer(o,delimiter=':')
		for m in range(len(bigAr)):
			params=bigAr[m]	
			fileWriter.writerow(params)
	else:
		t=csv.reader(open(infile,'rb'),delimiter=':')
		flines=[]
		ALL=[]
		for row in t:	
			try: f=1000*float(row[2])	#I had it in GHz.  Currently plotting in MHz.
			except: continue
			flines.append(1000*float(row[2]))
			ALL.append(row)
		flines=np.asarray(flines)
		if infile2 != 'None':
			tt= csv.reader(open(infile2,'rb'),delimiter=':')
			flines2=[]
			ALL2=[]
			for row in tt:
				try: f=1000*float(row[2])	#I had it in GHz.  Currently plotting in MHz.
				except: continue
				flines2.append(1000*float(row[2]))
				ALL2.append(row)
			flines2=np.asarray(flines2)
		if writeout == True or writeout == 'True':
			o=open(outfile,'wt')
			fileWriter=csv.writer(o,delimiter=':')
                        fileWriter.writerow(['Line #','Gaussian Center', 'crOffset from initial guess','Height','Width', 'Fit Measure','bin','Species0','Species', 'Line Rest Frequency', 'Transition','  ',' ','Species','Line Rest Frequency','Transition'])
			for m in range(len(bigAr)):
				dumar=[]
				params=bigAr[m]	
				dumar.extend(params)
				if infile =='recomb.csv':
					if abs(flines-1.5 - params[1]).min()<3:
						amin=np.argmin(abs(flines - params[1]))
						dumar.extend(ALL[amin])
					  	if plotQ==True or plotQ == 'True':
							try:
								pl.plot([params[1],params[1]],[-1*8,8],':',color='r')
							except:
								continue
					if infile2 != 'None':
						if abs(flines2-2. - params[1]).min()<2.5:
							amin=np.argmin(abs(flines2 - params[1]))
							dumar.extend(ALL2[amin])
					  		if plotQ==True or plotQ == 'True':
								pl.plot([params[1],params[1]],[-1*8,8],':',color='k')
				else: 
					if abs(flines - params[1]).min()<2.:
						amin=np.argmin(abs(flines - params[1]))
						dumar.extend(ALL[amin])
					  	if plotQ==True or plotQ == 'True':
							try:
								pl.plot([params[1],params[1]],[-1*8,8],':',color='r')
							except:
								continue
					if infile2 != 'None':
						if abs(flines2 - params[1]).min()<2.:
							amin=np.argmin(abs(flines2 - params[1]))
							dumar.extend(ALL2[amin])
					  		if plotQ== True or plotQ == 'True':
								pl.plot([params[1],params[1]],[-1*500,500],':',color='k')
				fileWriter.writerow(dumar)
		else:
			for m in range(len(bigAr)):
				dumar=[]
				params=bigAr[m]	
				dumar.extend(params)
				fGau = params[1]
				fOff = (64 - velcr)/c * fGau
				dF = delVel/c * fGau
				if abs(flines+fOff - fGau).min()<abs(dF):
					amin=np.argmin(abs(flines - fGau))
					dumar.extend(ALL[amin])
					fl =flines[amin]
				  	if plotQ== True or plotQ == 'True':
						pl.plot([fl,fl],[-1*500,500],':',color='r')
				if infile2 != 'None':
					if abs(flines2+fOff - fGau).min()<dF:
						amin=np.argmin(abs(flines2 - fGau))
						dumar.extend(ALL2[amin])
						fl =flines2[amin]
				  		if plotQ== True or plotQ == 'True':
							pl.plot([fl,fl],[-1*500,500],':',color='k')
	if writeout== True or writeout == 'True':			
		o.close()
