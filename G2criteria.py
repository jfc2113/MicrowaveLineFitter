### NOTE:  ANYWHERE YOU SEE (c-c0), YOU MAY NEED TO CHANGE THE CRITERIA FOR A DIFFERENT DATA SET.


def printFunc(quiet,text,text2='dum',text3='dum',text4='dum',text5='dum'):
	if quiet == False:
		print text
		if text2!='dum':
			print text2
			if text3!='dum':
				print text3
				if text4!='dum':
					print text4
					if text5!='dum':
						print text5

def c1(c,c0):
	if abs(c-c0)>0.3:
		return True
	else:
		return False 

def c2(w,BroadF):
	if abs(w)>BroadF:
		return True
	else:
		return False 

def c3(c,c0,rms,height,imrms):	
	if abs(height)/imrms >5 and height<0:
		if abs(c-c0) > 0.1 and abs(rms/height)>0.1 and rms/imrms > 1.:
	 		return True
		else:
			return False 
	else:
		if abs(c-c0) > 0.1 and abs(rms/height)>0.1 and rms/imrms > 1.1:
	 		return True
		return False 
def c4(rms,height):
	if abs(rms/height)>0.2:
		return True
	else: 	
		return False

def c5(c,c0,rms,imrms,PT_Height,q=False):
	if abs(PT_Height/imrms)>10:
	   if abs(rms)/imrms < 0.75:	##ADDED 26 AUG.  ...See if this fucks EVERYTHING up.
		return False
 	   if abs(c-c0)<0.6 and abs(rms/imrms) < .6: 
		return False
	   if abs(c-c0)<0.4 and abs(rms/PT_Height)< .10 and abs(rms/imrms) < 1:  #1.3 if I want to keep some of K6's beta & gamma recombs as 1comp, and some of L's recombs
		return False
	if abs(PT_Height/imrms)>7 and abs(PT_Height/imrms)<10:
	   if abs(rms)/imrms < 1.1:	##ADDED 26 AUG.
		return False
 	   if abs(c-c0)<0.6 and abs(rms/imrms) < 1:  	#.6 * freq/40000.
		return False
	   if abs(c-c0)<0.4 and abs(rms/PT_Height)< .12 and abs(rms/imrms) < 1.5:
		return False				##ADDED 26 AUG.
	if abs(PT_Height/imrms)<7:
	   if abs(rms)/imrms < 1.1:
		return False
	   if abs(c-c0)<0.6 and abs(rms/PT_Height)< .10 and abs(rms/imrms) < 1.5:   ##ADDED 26 AUG.
		return False

def c6(height,PT_Height,imrms,rms):
	if abs(PT_Height/imrms)> 10: 
  	   if abs(height-PT_Height) > imrms:
		if rms/imrms > .9:
		   return True
	if abs(PT_Height/imrms)> 7: 
  	   if abs(height-PT_Height) > imrms:
		if rms/imrms > 1:
		   return True

def c7(height,PT_Height,rms,imrms,q=False):
	printFunc(q,  "(height-PT_Height)/imrms:",(height-PT_Height)/imrms)
	printFunc(q,  "rms/imrms:",rms/imrms)
	printFunc(q,  "height/imrms:",height/imrms)
	if abs(height)/imrms > 1000: ### This is when a line is dynamic range limited.  Should only apply to maser lines in my data set.  Could reset based on your dynamic range limit.
		return False
	if abs(height-PT_Height)/imrms > 2 and rms/imrms > 2:# and abs(height/imrms) > 20:
		return True

	
	
def criteria(c,c0,w,rms,height,PT_Height,imrms,BroadF,q=False):

	if c1(c,c0)==True:
	   printFunc(q, "c1 == True")
	   if c6(height,PT_Height,imrms,rms) == True:
	  	printFunc(q,  "c6 == True")
		return True
	   else:
	    	if c5(c,c0,rms,imrms,PT_Height)==False:
	  	   printFunc(q, "c5 == False")
		else:
	  	   printFunc(q, "c5 == True")
		   return True
	if c2(w,BroadF)==True:
	   printFunc(q,  "c2 == True")
	   if c6(height,PT_Height,imrms,rms) == True:
	  	printFunc(q, "c6 == True")
		return True

	   else: 
		if c5(c,c0,rms,imrms,height)==False:
		  	printFunc(q, "c5 == False")
 	  	else:
	  		printFunc(q, "c5 == True")
			return True
	if c3(c,c0,rms,height,imrms)==True:
	   printFunc(q, "********* \n c3 == True \n ******")
  	   return True
	if c4(rms,height)==True:
	   printFunc(q, "In G2_criteria, c4 == True")
	   if c5(c,c0,rms,imrms,height)==False:
	  	printFunc(q, "In G2_criteria, c5 == False")
	   else:
	  	printFunc(q, "In G2_criteria, c5 == True")
		return True
	if c7(height,PT_Height,rms,imrms)==True:
	   print "In G2criteria, c7 == True"
	   if c5(c,c0,rms,imrms,height)==False:
	  	printFunc(q, "In G2_criteria, c5 == False")
	   else:
	  	printFunc(q, "In G2_criteria, c5 == True")
		return True
	return False
