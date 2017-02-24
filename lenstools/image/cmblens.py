"""

.. module:: cmblens
	:platform: Unix
	:synopsis: This module contains CMB lensing utilities


.. moduleauthor:: Andrea Petri <apetri@phys.columbia.edu>


"""

from abc import ABCMeta,abstractproperty,abstractmethod
import numpy as np
import astropy.units as u

try:
	import quicklens as ql
	ql = ql
except ImportError:
	ql = None

#####################
#Lens abstract class#
#####################

class Lens(object):

	__metaclass__ = ABCMeta

	#Make it a singleton
	_in_memory = None
	def __new__(cls,*args,**kwargs):
		if not(isinstance(cls._in_memory,cls)):
			cls._in_memory = object.__new__(cls,*args,**kwargs)
		return cls._in_memory

	##################
	#Abstract methods#
	##################

	@abstractmethod
	def __init__(self):
		pass

	@abstractmethod
	def generateTmap(angle,npixel):
		pass

	@abstractmethod
	def lensTmap(t,phifft,angle):
		pass

	@abstractmethod
	def phiTT(t,angle,powerTT,callback):
		pass

#########################
#quicklens functionality#
#########################

class QuickLens(Lens):

	#############
	#Constructor#
	#############
	
	def __init__(self):

		if ql is None:
			raise ImportError("This feature requires a quicklens installation!")

	######################################
	#Default CMB unlensed temperature map#
	######################################

	@staticmethod
	def generateTmap(angle,npixel,powerTT,callback):

		"""
		Default CMB temperature map in Fourier space

		"""

		#Calculate resolution
		resolution = angle.to(u.rad)/npixel

		#Parse the TT power spectrum
		if callback=="camb":
			
			if powerTT is None:
				Cl = ql.spec.get_camb_scalcl(lmax=3500)
			else:
				Cl = ql.spec.camb_clfile(powerTT)

		elif callback is None:
			raise NotImplementedError
		else:
			raise NotImplementedError

		#Build map
		pix = ql.maps.pix(npixel,resolution.value)
		teb_unl = ql.sims.tebfft(pix,Cl)

		#Return to user
		return teb_unl.tfft

	##################
	#Lens a CMB T map#
	##################

	@staticmethod
	def lensTmap(t,phifft,angle):
		
		"""
		Lens a CMB temperature map

		"""

		#Calculate number of pixels, resolution
		npixel = len(t)
		resolution = angle.to(u.rad)/npixel

		#Build TEB object, lens the map
		tqu_unl = ql.maps.tqumap(npixel,resolution.value,maps=[t,np.zeros_like(t),np.zeros_like(t)])
		tqu_len = ql.lens.make_lensed_map_flat_sky(tqu_unl,ql.maps.rfft(npixel,resolution.value,fft=phifft*resolution.value/npixel))

		#Return
		return tqu_len.tmap

	##################################
	#Quadratic TT potential estimator#
	##################################

	@staticmethod
	def phiTT(t,angle,powerTT,callback):

		"""
		Estimate the lensing potential with a quadratic TT estimator

		"""

		raise NotImplementedError

