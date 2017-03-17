"""

.. module:: cmblens
	:platform: Unix
	:synopsis: This module contains CMB lensing utilities


.. moduleauthor:: Andrea Petri <apetri@phys.columbia.edu>


"""

from abc import ABCMeta,abstractproperty,abstractmethod
import numpy as np
import astropy.units as u

from ..simulations.logs import logcmb

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
	_tcmb = 2.725*u.K

	#Make it a singleton
	_in_memory = None
	def __new__(cls,*args,**kwargs):
		if not(isinstance(cls._in_memory,cls)):
			cls._in_memory = object.__new__(cls,*args,**kwargs)
			cls._in_memory.resetCache()
		return cls._in_memory

	##############################################################

	#Set CMB temperature
	def setTCMB(self,t):
		assert(t.unit.physical_type==u"temperature")
		self._tcmb = t

	##############################################################

	###################
	#Noise in CMB maps#
	###################

	@staticmethod
	def _flat(ell,sigmaN):
		return (sigmaN**2)*np.ones_like(ell)

	@staticmethod
	def _detector(ell,sigmaN,fwhm,ellmax=None):
		power = (sigmaN**2)*np.exp(ell*(ell+1)*(fwhm**2)/(8*np.log(2)))
		if ellmax is not None:
			power[ell>ellmax] = 0.

		return power

	##############################################################

	#Reset the cache
	def resetCache(self):
		self._cache = dict()
		self._cache["angle"] = -1
		self._cache["npixel"] = -1
		self._cache["lmax"] = -1

	#Build multipoles cache
	def buildEllCache(self,angle,npixel,lmax):

		if (angle!=self._cache["angle"]) or (npixel!=self._cache["npixel"]) or (lmax!=self._cache["lmax"]):

			#Log
			logcmb.debug("Building multipole cache...")
			
			#Angle
			self._cache["angle"] = angle
			self._cache["npixel"] = npixel
			self._cache["lmax"] = lmax

			#Multipoles
			f = np.fft.fftfreq(npixel)			
			ell_x,ell_y = np.meshgrid(f,f,indexing="ij")
			self._cache["ell2"] = (ell_x**2 + ell_y**2)*((2.0*np.pi*npixel/(angle.to(u.rad).value))**2)
			self._cache["ell2"][0,0] = 1.0
			self._cache["ell"] = np.sqrt(self._cache["ell2"])

	#Build TT unlensed power spectrum cache
	def buildUnlTTCache(self,powerTT,callback):

		#Log
		logcmb.debug("Building unlensed CMB power spectrum cache...")

		#Cache name
		self._cache["powerTT_unl_name"] = str(powerTT)

		#Load power spectra
		self._cache["Cl_unl"] = self.getPower(powerTT,callback,self._cache["lmax"])
		self._cache["powerTT_unl"] = self.extractTT(self._cache["Cl_unl"])

	#Build TT lensed power spectrum cache
	def buildLensedTTCache(self,powerTT,callback,noise_keys):

		#Log
		logcmb.debug("Building lensed CMB power spectrum cache...")

		#Cache name
		self._cache["powerTT_lensed_name"] = str(powerTT)

		#Load power spectra
		self._cache["Cl_lensed"] = self.getPower(powerTT,callback,self._cache["lmax"])
		ClTT_lensed = self.extractTT(self._cache["Cl_lensed"]) 
		self._cache["powerTT_lensed"] = ClTT_lensed

		#Build inverse variance filter: add noise component to the observed power spectrum
		self._cache["powerTT_obs"] = ClTT_lensed.copy()
		if noise_keys is not None:
			self._cache["powerTT_obs"] += self.getNoise(np.arange(len(ClTT_lensed)),noise_keys)
		
		#Construct inverse variance filter
		self._cache["1/powerTT_obs"] = 1./self._cache["powerTT_obs"]
		self._cache["1/powerTT_obs"][[0,1]] = 0.

	#Noise
	def getNoise(self,ell,noise_keys):
		assert "kind" in noise_keys,"Format of the noise keys must be {'kind':'white,detector','sigmaN':value,'fwhm':value}"

		if noise_keys["kind"]=="white":

			#White noise
			sigmaN = noise_keys["sigmaN"]
			return self._flat(ell,sigmaN.to(u.uK*u.rad).value)

		elif noise_keys["kind"]=="detector":

			#Detector noise
			sigmaN = noise_keys["sigmaN"]
			fwhm = noise_keys["fwhm"]
			return self._detector(ell,sigmaN.to(u.uK*u.rad).value,fwhm.to(u.rad).value)

		else:
			raise NotImplementedError("Noise kind '{0}' not implemented: choose (white/detector)".format(noise_keys["kind"]))

	##############################################################

	##################
	#Abstract methods#
	##################

	@abstractmethod
	def __init__(self):
		pass

	@abstractmethod
	def getPower(self,power,callback,lmax):
		pass

	@abstractmethod
	def extractTT(Cl):
		pass

	@abstractmethod
	def buildNormCache(self,npixel,resolution,estimator,Tfilter):
		pass

	@abstractmethod
	def buildQQCache(self,npixel,resolution,estimator,Tfilter):
		pass

	@abstractmethod
	def generateTmap(angle,npixel):
		pass

	@abstractmethod
	def lensTmap(self,t,angleT,kappa,angleKappa):
		pass

	@abstractmethod
	def phiTT(tfft,angle,powerTT,callback,noise_keys,lmax,filtering):
		pass

	@abstractmethod
	def N0TT(self,l_edges,angle,npixel,powerTT,callback,noise_keys,lmax):
		pass

#########################
#quicklens functionality#
#########################

class QuickLens(Lens):

	#############
	#Constructor#
	#############
	
	def __init__(self):

		#Quicklens
		if ql is None:
			raise ImportError("This feature requires a quicklens installation!")

	def getPower(self,power,callback,lmax):
		
		#Select correct format
		if callback in ("camb_dimensionless","camb_uk"):
			
			if power is None:
				Cl = ql.spec.get_camb_scalcl(lmax=3500)
			else:
				Cl = ql.spec.camb_clfile(power,lmax=lmax)

			#Scale to units to uK^2 if loaded power spectrum is dimensionless
			if callback=="camb_dimensionless":
				tcmb = self._tcmb.to(u.uK).value

				for s in ("cltt","clee","clbb","clte","clpp","cltp"):
					if hasattr(Cl,s):
						logcmb.debug("Scaling loaded {0} coefficients to uK^2 units...".format(s))
						cl = getattr(Cl,s)
						cl *= tcmb**2

		elif callback is None:
			raise NotImplementedError
		else:
			raise NotImplementedError

		#Return
		return Cl

	@staticmethod
	def extractTT(Cl):
		return Cl.cltt

	#Build normalization cache
	def buildNormCache(self,npixel,resolution,estimator,Tfilter):

		#Log
		logcmb.debug("Building quadratic estimator normalization cache...")

		self._cache["norm"] = estimator.fill_resp(estimator,ql.maps.cfft(npixel,resolution),Tfilter,Tfilter)
		self._cache["norm"].fft[0,0] = 1.0

		#Log
		logcmb.debug("Normalization cache complete")

	#Build QQ cache
	def buildQQCache(self,npixel,resolution,estimator,Tfilter):
		self._cache["qq"] = estimator.fill_clqq(ql.maps.cfft(npixel,resolution),Tfilter,Tfilter,Tfilter)

	######################################
	#Default CMB unlensed temperature map#
	######################################

	def generateTmap(self,angle,npixel,powerTT,callback,lmax):

		"""
		Default CMB temperature map in Fourier space

		"""

		#Calculate resolution
		resolution = angle.to(u.rad)/npixel

		#Parse the TT power spectrum
		if ("powerTT_unl" not in self._cache) or (powerTT!=self._cache["powerTT_unl_name"]):
			self.buildEllCache(angle,npixel,lmax)
			self.buildUnlTTCache(powerTT,callback)

		#Build map
		pix = ql.maps.pix(npixel,resolution.value)
		teb_unl = ql.sims.tebfft(pix,self._cache["Cl_unl"])

		#Return to user
		return teb_unl.tfft*npixel/resolution.value

	##################
	#Lens a CMB T map#
	##################

	def lensTmap(self,t,angleT,kappa,angleKappa):
		
		"""
		Lens a CMB temperature map

		"""

		#Calculate number of pixels, resolution
		npixelT = len(t)
		resolutionT = angleT.to(u.rad)/npixelT
		npixelPhi = len(kappa)
		resolutionPhi = angleKappa.to(u.rad)/npixelPhi

		#Build TEB object
		tqu_unl = ql.maps.tqumap(npixelT,resolutionT.value,maps=[t,np.zeros_like(t),np.zeros_like(t)])
		
		#Compute potential
		kappa_fft = ql.maps.rmap(npixelPhi,resolutionPhi.value,map=kappa).get_rfft()
		phi_fft = ql.maps.rfft(npixelPhi,resolutionPhi.value,fft=kappa_fft.fft * 2.0 / self._cache["ell2"][:,:npixelPhi//2+1])

		#Zero out high multipoles
		phi_fft.fft[self._cache["ell"][:,:npixelPhi//2+1]>self._cache["lmax"]] = 0.

		#Lens the map
		tqu_len = ql.lens.make_lensed_map_flat_sky(tqu_unl,phi_fft)

		#Return
		return tqu_len.tmap

	##################################
	#Quadratic TT potential estimator#
	##################################

	def phiTT(self,tfft,angle,powerTT,callback,noise_keys,lmax,filtering):

		"""
		Estimate the lensing potential with a quadratic TT estimator

		"""

		#Calculate number of pixels, resolution
		npixel = len(tfft)
		resolution = angle.to(u.rad)/npixel

		#######################################################
		#Build the caches for ell,C_ell if not present already#
		#######################################################

		if ("powerTT_obs" not in self._cache) or (str(powerTT)!=self._cache["powerTT_lensed_name"]) or (lmax!=self._cache["lmax"]):
			self.buildEllCache(angle,npixel,lmax)
			self.buildLensedTTCache(powerTT,callback,noise_keys)
			recompute_norm = True
		else:
			recompute_norm = False
		
		#Build the TT estimator
		estimator = ql.qest.lens.phi_TT(self._cache["powerTT_lensed"])

		#Apply inverse variance filter to the observation
		tfft = ql.maps.cfft(npixel,resolution.value,fft=tfft*resolution.value/npixel)
		tfft = tfft * self._cache["1/powerTT_obs"]

		#Evaluate the estimator on kappa and its normalization
		phi_eval = estimator.eval_flatsky(tfft,tfft,npad=1)

		if recompute_norm:
			self.buildNormCache(npixel,resolution.value,estimator,self._cache["1/powerTT_obs"])

		#Normalize the estimator
		phifft = phi_eval/self._cache["norm"]

		#Apply filter
		if filtering is not None:

			if filtering=="wiener":
				phifft = phifft * (self._cache["powerTT_lensed"] * self._cache["1/powerTT_obs"])
			else:
				phifft.fft *= filtering(self._cache["ell"])

		#Return
		return phifft.fft*npixel/resolution.value

	##################
	#Analytic N0 bias#
	##################

	def N0TT(self,l_edges,angle,npixel,powerTT,callback,noise_keys,lmax): 

		#######################################################
		#Build the caches for ell,C_ell if not present already#
		#######################################################

		if ("powerTT_obs" not in self._cache) or (str(powerTT)!=self._cache["powerTT_lensed_name"]) or (lmax!=self._cache["lmax"]):
			self.buildEllCache(angle,npixel,lmax)
			self.buildLensedTTCache(powerTT,callback,noise_keys)
			recompute_norm = True
			recompute_qq = True
		else:
			recompute_norm = False
			recompute_qq = True

		#Get a handle on the estimator
		resolution = angle.to(u.rad).value/npixel
		estimator = ql.qest.lens.phi_TT(self._cache["powerTT_lensed"])

		#Recompute normalization
		if ("norm" not in self._cache) or recompute_norm:
			self.buildNormCache(npixel,resolution,estimator,self._cache["1/powerTT_obs"])

		#Recompute qq
		if ("qq" not in self._cache) or recompute_qq:
			self.buildQQCache(npixel,resolution,estimator,self._cache["1/powerTT_obs"])

		#Get normalized QQ estimate
		nlqq = self._cache["qq"] / self._cache["norm"]**2

		#Azimuthal average
		nlqq_azimuthal = nlqq.get_ml(l_edges)

		#Return
		return nlqq_azimuthal.ls,nlqq_azimuthal.specs["cl"].real


