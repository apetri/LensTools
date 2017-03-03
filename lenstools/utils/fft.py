from __future__ import division

from abc import ABCMeta,abstractproperty,abstractmethod

import numpy as np

##############################################
###########FFTEngine abstract class###########
##############################################

class FFTEngine(object):

	__metaclass__ = ABCMeta

	"""
	Class handler of Fourier transforms needed for lenstools computations

	"""

	#####################################################################
	######################Abstract methods###############################
	#####################################################################

	@abstractmethod
	def fft2(self,x):
		pass

	@abstractmethod
	def ifft2(self,x):
		pass

	@abstractmethod
	def rfft2(self,x):
		pass

	@abstractmethod
	def irfft2(self,x):
		pass

	@abstractmethod
	def rfftn(self,x):
		pass

	@abstractmethod
	def irfftn(self,x):
		pass

	###################################################################################
	######################Default, non--abstract methods###############################
	###################################################################################

	def fftfreq(self,n):
		return np.fft.fftfreq(n)

	def rfftfreq(self,n,d=1.0):
		
		if not (isinstance(n,int)):
			raise ValueError("n should be an integer")
    	
		val = 1.0/(n*d)
		n_half = n//2 + 1
		results = np.arange(0, n_half, dtype=int)
		return results * val
		


##############################################
###########NUMPYFFTPack class#################
##############################################

class NUMPYFFTPack(FFTEngine):

	###############################################################################################
	#########################Abstract methods implementation#######################################
	###############################################################################################

	def fft2(self,x):
		return np.fft.fft2(x)

	def ifft2(self,x):
		return np.fft.ifft2(x)

	def rfft2(self,x):
		return np.fft.rfft2(x)

	def irfft2(self,x):
		return np.fft.irfft2(x)

	def rfftn(self,x):
		return np.fft.rfftn(x)

	def irfftn(self,x):
		return np.fft.irfftn(x)
