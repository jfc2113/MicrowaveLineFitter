import numpy.fft as fft, numpy as np, pylab as pl


def padding(smallFreq,smallFlux,plotQ='False'):
	fftFlux=fft.fftshift(fft.fft(smallFlux))

	tb=np.zeros((len(fftFlux)))

	fftFlux=np.append(tb,fftFlux,axis=0)
	fftFlux=np.append(fftFlux,tb,axis=0)
	
	InterpFlux=fft.ifft(fft.ifftshift(fftFlux))
	InterpFlux=InterpFlux[0:len(InterpFlux)-2]
	#InterpFreq=np.linspace(smallFreq[0],smallFreq[len(smallFreq)-1],len(InterpFlux))-1 *(abs(smallFreq[0]-smallFreq[len(smallFreq)-1]))/len(InterpFlux)
	InterpFreq=np.linspace(smallFreq[0],smallFreq[len(smallFreq)-1],len(InterpFlux))
	if plotQ == 'True':
		pl.plot(smallFreq,smallFlux,label='original')
		pl.plot(InterpFreq,3*np.real(InterpFlux),label='interp')
		pl.legend()
		pl.show()
	return InterpFreq,3*np.real(InterpFlux)




###This is shifted correctly.  It's not normalized to the same total flux though.  
