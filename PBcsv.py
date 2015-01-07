import numpy as np, pylab as pl, csv

def PBcor(reg,fAR,fdAR):
	if reg =='K6':
		gg,ee,dd,cc,bb,aa = [  8.41602173e-24,  -1.53590760e-18,   1.10043699e-13,  -3.87260572e-09,   6.71852727e-05,   5.41181354e-01]
	if reg == 'L':
		gg,ee,dd,cc,bb,aa = [ -2.17496751e-21,   4.29989569e-16,  -3.35563412e-11,   1.29344001e-06,  -2.46042047e-02,   1.85986566e+02]
	if reg =='K4':
		gg,ee,dd,cc,bb,aa = [ -2.55236267e-23,   5.08709478e-18,  -4.02860315e-13,   1.58227948e-08,  -3.06457861e-04,   3.34893004e+00]
	PBrat = aa + bb*fAR + cc*fAR**2 + dd*fAR**3 + ee*fAR**4 + gg*fAR**5
	return fdAR * PBrat

def PBcor_oneLine(reg,f,h):
	if reg =='K6':
		gg,ee,dd,cc,bb,aa = [  8.41602173e-24,  -1.53590760e-18,   1.10043699e-13,  -3.87260572e-09,   6.71852727e-05,   5.41181354e-01]
	if reg == 'L':
		gg,ee,dd,cc,bb,aa = [ -2.17496751e-21,   4.29989569e-16,  -3.35563412e-11,   1.29344001e-06,  -2.46042047e-02,   1.85986566e+02]
	if reg =='K4':
		gg,ee,dd,cc,bb,aa = [ -2.55236267e-23,   5.08709478e-18,  -4.02860315e-13,   1.58227948e-08,  -3.06457861e-04,   3.34893004e+00]
	PBrat = aa + bb*f + cc*f**2 + dd*f**3 + ee*f**4 + gg*f**5
	return h * PBrat


def CSVpb(reg,infile,outfile):
	t=csv.reader(open(infile,'rt'),delimiter=':')
	o=open(outfile,'wt')
	fileWriter=csv.writer(o,delimiter=':')
	fileWriter.writerow(['Line Rest Frequency (MHz)', 'Species', 'Transition', 'Velocity (km/s)','Height','Width (km/s)','Line Center (MHz)','Line Type'])
	for row in t:	
		print "row= first",row
		print "row[0] =", row[0]
		try: f=float(row[6])	
		except: continue
		print "row= before",row
		h = float(row[4])
		row[4] = PBcor_oneLine(reg,f,h)
		print "row after",row
		fileWriter.writerow(row)
	o.close()
