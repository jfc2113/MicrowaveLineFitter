import numpy as np, pylab as pl
from math import *



cutchans=5		##This might work
print "Importing LineFind.  Cutchans = ", cutchans

def sd_MAD(yar):
	med = np.median(yar)
	return 1.4826*np.median(np.abs(yar-med))

def cut(xar,yar,cutSig,Set_rms = 'dummy'):
	fullcutX = np.zeros((2))
	fullcutY = np.zeros((2))
	rmss=[]
	iii=0
	for freqCenter in np.arange(30750.,50000.,1850.):
		xmask = abs(xar-freqCenter)<925.
		xnow = xar[xmask]
		ynow = yar[xmask]
		if len(xnow) ==0:
			continue
		if Set_rms =='dummy':
			rms = sd_MAD(ynow)
		else: rms = Set_rms[iii]
		mycut=abs(ynow)<cutSig*rms
		fullcutX=np.concatenate((fullcutX,xnow[mycut]))
		fullcutY=np.concatenate((fullcutY,ynow[mycut]))
		rmss.append(rms)
		iii+=1
	print "rms array!",rmss
	return fullcutX[2:],fullcutY[2:],rmss

def backwardsCut(mainii,cc,xar,xar_old,del_list=[]):
	for dumindx in np.arange(-1*cc,0,1):
		try: item2 = xar_old[mainii+dumindx]
		except: continue		#This should kill any channels at the edges of the spectrum.
		if item2 in xar:
			jjj=np.argmin(abs(xar-item2))
			del_list.append(jjj)
	return del_list
def forwardsCut(mainii,cc,xar,xar_old,del_list=[]):
	pushforward=0
	for dumindx in np.arange(0,cc+1,1):
		try: item2 = xar_old[mainii+dumindx]
		except: continue		#This should kill any channels at the edges of the spectrum.
		if item2 in xar:
			jjj=np.argmin(abs(xar-item2))
			del_list.append(jjj)
		else: 
			pushforward=dumindx
	return pushforward, del_list

		
def neighborCut(xar,xar_old,yar,yar_old,cc): 	#cc is the number of channels you want to cut on each side.
	mainii=0
	dumVar=0
	while mainii<len(xar_old):
		item = xar_old[mainii]
		if not item in xar:	
			DL=backwardsCut(mainii,cc,xar,xar_old,del_list=[])
			PF,DL = forwardsCut(mainii,cc,xar,xar_old,del_list=DL)
			while PF>0:	
				mainii+=PF	
				PF,DL = forwardsCut(mainii,cc,xar,xar_old,del_list=DL)
			xar=np.delete(xar,DL)
			yar=np.delete(yar,DL)
			mainii+=cc
			#print " "
		mainii+=1	
	return xar, yar

def spect(infile):
	t=open(infile)
	freq=[]
	fluxD=[]

	for line in t:
		l=line.split()
		try: a=float(l[0])
		except:	continue
		b=float(l[1])
		freq.append(a)
		fluxD.append(b)
	t.close()
	return np.asarray(freq),np.asarray(fluxD)



def LineFind(f,fd,r_rms='dummy',plotQ='False',cc1=cutchans,rmsthresh0=4,rmsthresh1=3.0):
	fLF=f	#frequency, line free
	fdLF=fd	#flux dens, line free
	fLF,fdLF,rmsfinal=cut(fLF,fdLF,rmsthresh0,r_rms)		##Do this to evaluate the standard deviation using the Median Absolute Deviation.
	fLF,fdLF=neighborCut(fLF,f,fdLF,fd,2)		##Cut 2 channels before evaluating the standard deviation again.
	fLF,fdLF,rmsff=cut(fLF,fdLF,rmsthresh1,r_rms)		##Evaluate the standard deviation using only the non-line data.
	fLF,fdLF,rmsfinal=cut(f,fd,rmsthresh1,rmsff)		##Apply the cut to the original data using the calculated rms from the non-line data
	fLF,fdLF=neighborCut(fLF,f,fdLF,fd,cc1)	##Neighbor cut!
 	freq_diff = []
 	fd_diff = []
	for ii in range(len(f)):
		item=f[ii]
		if not item in fLF:
			freq_diff.append(item)
			fd_diff.append(fd[ii])
	if r_rms!='dummy':
		rmsfinal=r_rms
	#print len(f),len(fd)
	#print len(freq_diff),len(fd_diff)
	if plotQ=='True':
		fig = pl.figure()
		ax = fig.add_subplot(111)
		pl.scatter(f,fd,s=0.5,color='c',marker='+')
		pl.scatter(freq_diff,fd_diff,s=1.0,marker='o',color='y')
		pl.legend()
		pl.show()
	return np.asarray(freq_diff),np.asarray(fd_diff),rmsfinal
