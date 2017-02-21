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

	@abstractmethod
	def __init__(self):
		pass

	@abstractmethod
	def defaultTmap(angle,npixel):
		pass

	@abstractmethod
	def lensTmap(t,phifft,npixel,angle):
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

	##########################
	#Default CMB unlensed map#
	##########################

	@staticmethod
	def defaultTmap(angle,npixel):

		"""
		Default CMB temperature map in Fourier space

		"""

		#Calculate resolution
		resolution = angle.to(u.rad)/npixel
		lmax = 3000

		#Build map
		Cl = ql.spec.get_camb_scalcl(lmax=lmax)
		pix = ql.maps.pix(npixel,resolution.value)
		teb_unl = ql.sims.tebfft(pix,Cl)

		#Return to user
		return teb_unl.tfft

	################
	#Lens a CMB map#
	################

	@staticmethod
	def lensTmap(t,phifft,npixel,angle):
		
		"""
		Lens a CMB temperature map

		"""

		#Calculate resolution
		resolution = angle.to(u.rad)/npixel

		#Build TEB object, lens the map
		tqu_unl = ql.maps.tqumap(npixel,resolution.value,maps=[t,np.zeros_like(t),np.zeros_like(t)])
		tqu_len = ql.lens.make_lensed_map_flat_sky(tqu_unl,ql.maps.rfft(npixel,resolution.value,fft=phifft))

		#Return
		return tqu_len.tmap

