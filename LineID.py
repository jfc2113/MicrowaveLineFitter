import numpy as np, pylab as pl, csv
from math import *

###  This compares the data in bigAr, which is based on the output of Gfit.Gaufit and managed through FullCode.py, to data in one or two csv files containing line information.  
### If writeout == True (or 'True'), it writes a csv file 

def lineCompareAll(outfile,bigAr,infile1='None',infile2 = 'None',plotQ=False,writeout=True,inFileUnits = 'GHz',velshift=64,velRecomb=75, dV_recomb=20, velMol = 70, dV_mol=17,best=-9999,plotlineHeight='dum'):
	if plotlineHeight=='dum':
		d0,d1,d2,heights,d4,d5,binCodes = zip(*bigAr)
		mask=np.asarray(binCodes)<1.5
		heights=np.absolute(np.asarray(heights)[mask])
		heights = heights[np.argsort(heights)]
		plotlineHeight=heights.max()
	if best==-9999:
		best=velshift
	c = 2.9979E5 #  in km/s	
	if infile2 =='recomb.csv':
		infile2=infile
		infile = 'recomb.csv'
	flines=[]

	### Fill arrays containing the line data from the csv files. 
	if infile1 != 'None' and infile1 !='dum' and  infile1 != 'none' and infile1 !='Dum':
		t=csv.reader(open(infile1,'rb'),delimiter=':')
		ALL=[]
		for row in t:	
			try: f=1000*float(row[2])	
			except: continue
			if inFileUnits == 'GHz':
				flines.append(1000*float(row[2]))
			else:
				flines.append(float(row[2]))
			ALL.append(row)
		flines=np.asarray(flines)
		if len(flines)==0:
			print "Seems like infile1 is incorrectly formatted.  Can't proceed properly.  Recommend reformatting to a text csv file with field delimiter = ':'."
	if infile2 != 'None' and infile2 !='dum':
		tt= csv.reader(open(infile2,'rb'),delimiter=':')
		flines2=[]
		ALL2=[]
		for row in tt:
			try: f=1000*float(row[2])	#I had it in GHz.  Currently plotting in MHz.
			except: continue
			if inFileUnits == 'GHz':
				flines2.append(1000*float(row[2]))
			else:
				flines2.append(float(row[2]))
			ALL2.append(row)
		flines2=np.asarray(flines2)
		if len(flines2)==0:
			print "Seems like infile2 is incorrectly formatted.  Can't proceed properly.  Recommend reformatting. to a text csv file with field delimiter = ':'."

	## Compare the line data to the detected line frequencies to determine line-identifications.
	if writeout == True or writeout == 'True':
		## Make a csv file output.
		o=open(outfile,'wt')
		fileWriter=csv.writer(o,delimiter=':')
                fileWriter.writerow(['Line #','Gaussian Center', 'crOffset from initial guess','Height','Width', 'Fit Measure','bin','Species','Molecule', 'Line Rest Frequency', 'Transition','Line Strength','Line List','Species','Molecule','Line Rest Frequency (%s)' % inFileUnits,'Transition','Line Strength','Line List'])
	for m in range(len(bigAr)):		#bigAr contains the detected gaussian fit parameters.
		dumar=[]
		params=bigAr[m]	
		dumar.extend(params)
		if len(flines)==0:
			if writeout == True or writeout == 'True':
				fileWriter.writerow(dumar)
				continue
		if infile1 =='recomb.csv' or ALL[0][4]=='Recomb' or ALL[0][5]=='Recomb' or ALL[0][4]=='recomb' or ALL[0][5]=='recomb':
			vlines = velshift + (flines-params[1])/params[1] * c
			if abs(vlines-velRecomb).min() < dV_recomb:
				nowmask = abs(vlines-velRecomb) < dV_recomb
				vlinesnow=vlines[nowmask]
				for i in range(len(vlinesnow)):
					aimin = np.argmin(abs(vlinesnow-velRecomb))
					amin = np.argmin(abs(vlines-vlinesnow[aimin]))
					dumar.extend(ALL[amin][:4])
					dumar.extend(['--'])
					if ALL[0][4]=='Recomb' or ALL[0][4]=='recomb':
						dumar.extend(ALL[amin][4:5])
					else:
						if ALL[0][5]=='Recomb' or ALL[0][5:6]=='recomb':
							dumar.extend(ALL[amin][5])
					fl = flines[amin]
					vlinesnow = np.concatenate((vlinesnow[:aimin],vlinesnow[aimin+1:]))
				  	if plotQ==True or plotQ == 'True':
						flp= fl*(1+(velshift-velRecomb)/c)
						try:
							pl.plot([flp,flp],[-1*plotlineHeight,plotlineHeight],':',color='r')
						except:
							continue
		else:
			vlines = velshift  +  (flines-params[1])/params[1] * c
			if abs(vlines-velMol).min() < dV_mol:
				nowmask = abs(vlines-velMol) < dV_mol
				vlinesnow=vlines[nowmask]
				for i in range(len(vlinesnow)):
					aimin = np.argmin(abs(vlinesnow-best))
					amin = np.argmin(abs(vlines-vlinesnow[aimin]))
					dumar.extend(ALL[amin][:6])
					fl = flines[amin]
					vlinesnow = np.concatenate((vlinesnow[:aimin],vlinesnow[aimin+1:]))
				  	if plotQ==True or plotQ == 'True':
						flp= fl*(1+(velshift-best)/c)
							
						try:
							pl.plot([flp,flp],[-1*plotlineHeight,plotlineHeight],':',color='r')
						except:
							continue
			
		if infile2 != 'None' and infile2 !='dum':
			if len(flines2)==0:
				continue
			vlines2 = velshift + (flines2-params[1])/params[1] * c
			if abs(vlines2-velMol).min() < dV_mol:
				nowmask = abs(vlines2-velMol) < dV_mol
				vlinesnow=vlines2[nowmask]
				for i in range(len(vlinesnow)):
					aimin = np.argmin(abs(vlinesnow-best))
					amin = np.argmin(abs(vlines2-vlinesnow[aimin]))
					dumar.extend(ALL2[amin][:6])
					fl = flines2[amin]
					vlinesnow = np.concatenate((vlinesnow[:aimin],vlinesnow[aimin+1:]))
				  	if plotQ==True or plotQ == 'True':
						try:
							flp= fl*(1+(velshift-best)/c)
							pl.plot([flp,flp],[-1*plotlineHeight,plotlineHeight],':',color='k')
						except:
							continue
		if writeout == True or writeout == 'True':
			fileWriter.writerow(dumar)

	if writeout== True or writeout == 'True':			
		o.close()
