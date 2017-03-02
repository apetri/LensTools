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
			cls._in_memory.reset()
		return cls._in_memory

	##############################################################

	#Reset the cache
	def reset(self):
		self._cache = dict()
		self._cache["angle"] = -1
		self._cache["npixel"] = -1

	#Build multipoles cache
	def buildEllCache(self,angle,npixel):

		if (angle!=self._cache["angle"]) or (npixel!=self._cache["npixel"]):
			
			#Angle
			self._cache["angle"] = angle
			self._cache["npixel"] = npixel

			#Multipoles			
			ell_x,ell_y = np.meshgrid(np.fft.fftfreq(npixel),np.fft.fftfreq(npixel),indexing="ij")
			self._cache["ell2"] = (ell_x**2 + ell_y**2)*((2.0*np.pi*npixel/(angle.to(u.rad).value))**2)
			self._cache["ell"] = np.sqrt(self._cache["ell2"])

	##############################################################

	##################
	#Abstract methods#
	##################

	@abstractmethod
	def __init__(self):
		pass

	@abstractmethod
	def getPower(power,callback):
		pass

	@abstractmethod
	def buildPowerTTCache(self,powerTT_th,powerTT_obs,callback):
		pass

	@abstractmethod
	def generateTmap(angle,npixel):
		pass

	@abstractmethod
	def lensTmap(t,phifft,angle):
		pass

	@abstractmethod
	def phiTT(tfft,angle,powerTT_th,powerTT_obs,callback):
		pass

#########################
#quicklens functionality#
#########################

class QuickLens(Lens):

	_tiny = 1.0e-24

	#############
	#Constructor#
	#############
	
	def __init__(self):

		#Quicklens
		if ql is None:
			raise ImportError("This feature requires a quicklens installation!")

	@staticmethod
	def getPower(power,callback):
		
		#Select correct format
		if callback=="camb":
			
			if power is None:
				Cl = ql.spec.get_camb_scalcl(lmax=3500)
			else:
				Cl = ql.spec.camb_clfile(power)

		elif callback is None:
			raise NotImplementedError
		else:
			raise NotImplementedError

		#Return
		return Cl

	#Build TT power spectrum cache
	def buildPowerTTCache(self,powerTT_th,powerTT_obs,callback):

		#Npixel, resolution
		npixel = self._cache["npixel"]
		resolution = self._cache["angle"].to(u.rad).value/npixel

		#Load power spectra
		self._cache["ClTT_th"] = self.getPower(powerTT_th,callback).cltt
		ClTT_th = self._cache["ClTT_th"] 
		
		if powerTT_obs is None:
			ClTT_obs = None
		else:
			ClTT_obs = self.getPower(powerTT_obs,callback).cltt 

		#Build the power TT caches
		self._cache["powerTT_th"] = np.interp(self._cache["ell"].flatten(),np.arange(len(ClTT_th)),ClTT_th).reshape(self._cache["ell"].shape)

		if ClTT_obs is None:
			self._cache["powerTT_obs"] = self._cache["powerTT_th"]
		else:
			self._cache["powerTT_obs"] = np.interp(self._cache["ell"].flatten(),np.arange(len(ClTT_obs)),ClTT_obs).reshape(self._cache["ell"].shape)

		self._cache["powerTT_obs"][0,0] = self._tiny
		self._cache["1/powerTT_obs"] = 1./self._cache["powerTT_obs"]
		self._cache["1/powerTT_obs"][0,0] = self._tiny

		#Build the cache for the normalization factor
		self._cache["norm"] = ql.maps.cfft(npixel,resolution,fft=np.zeros((npixel,npixel),dtype=np.complex))

	######################################
	#Default CMB unlensed temperature map#
	######################################

	def generateTmap(self,angle,npixel,powerTT,callback):

		"""
		Default CMB temperature map in Fourier space

		"""

		#Calculate resolution
		resolution = angle.to(u.rad)/npixel

		#############################################################
		#TODO: Build the caches for ell,C_ell if not present already#
		#############################################################

		#Parse the TT power spectrum
		Cl = self.getPower(powerTT,callback)

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

		#############################################################
		#TODO: Build the caches for ell,C_ell if not present already#
		#############################################################

		#Build TEB object, lens the map
		tqu_unl = ql.maps.tqumap(npixel,resolution.value,maps=[t,np.zeros_like(t),np.zeros_like(t)])
		tqu_len = ql.lens.make_lensed_map_flat_sky(tqu_unl,ql.maps.rfft(npixel,resolution.value,fft=phifft*resolution.value/npixel))

		#Return
		return tqu_len.tmap

	##################################
	#Quadratic TT potential estimator#
	##################################

	def phiTT(self,tfft,angle,powerTT_th,powerTT_obs,callback):

		"""
		Estimate the lensing potential with a quadratic TT estimator

		"""

		#Calculate number of pixels, resolution
		npixel = len(tfft)
		resolution = angle.to(u.rad)/npixel

		#######################################################
		#Build the caches for ell,C_ell if not present already#
		#######################################################

		if "powerTT_th" not in self._cache:
			self.buildEllCache(angle,npixel)
			self.buildPowerTTCache(powerTT_th,powerTT_obs,callback)
		
		#Build the TT estimator
		estimator = ql.qest.lens.phi_TT(self._cache["ClTT_th"])

		#Apply inverse variance filter to the observation
		tfft = tfft * self._cache["1/powerTT_obs"][:,:npixel//2+1] * (resolution.value**2)

		#Evaluate the estimator on kappa and its normalization
		phi_eval = estimator.eval(ql.maps.rfft(npixel,resolution.value,fft=tfft))
		phi_norm = estimator.fill_resp(estimator,self._cache["norm"],self._cache["1/powerTT_obs"],self._cache["1/powerTT_obs"])
		phi_norm.fft[0,0] = 1.

		#Return the normalized estimator
		phifft = phi_eval/phi_norm
		return phifft.fft 

