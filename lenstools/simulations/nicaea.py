import numpy as np
from astropy.cosmology import w0waCDM

try:
	from ..extern import _nicaea
	_nicaea=_nicaea
except ImportError:
	_nicaea=None


##########################################
##########NicaeaSettings class############
##########################################

class NicaeaSettings(dict):

	"""
	Class handler of the code settings (non linear modeling, tomography, transfer function, etc...)
	"""

	def __init__(self):

		super(NicaeaSettings,self).__init__()
		
		#Default settings
		self["snonlinear"]="smith03"
		self["stransfer"]="eisenhu"
		self["sgrowth"]="growth_de"
		self["sde_param"]="linder"
		self["normmode"]="norm_s8"
		self["stomo"]="tomo_all"
		self["sreduced"]="none"
		self["q_mag_size"]=1.0

	@classmethod
	def default(cls):

		"""
		Generate default settings

		:returns: NicaeaSettings defaults instance

		"""

		return cls()

	@property
	def knobs(self):

		"""
		Lists available settings to tune

		"""

		return self.keys()

	def available(self,knob):

		"""
		Given a settings, lists all the possible values
		
		"""
		
		#Available settings
		if knob=="snonlinear":
			return [ "linear", "pd96", "smith03", "smith03_de", "coyote10", "halodm" ]
		elif knob=="stransfer":
			return [ "bbks", "eisenhu", "eisenhu_osc", "camb_vinschter", "camb", "be84" ]
		elif knob=="sgrowth":
			return [ "heath", "growth_de", "camb_vinschter_gr" ]
		elif knob=="sde_param":
			return [ "jassal", "linder", "earlyDE", "poly_DE" ]
		elif knob=="normmode":
			return [ "norm_s8" , "norm_as" ]
		elif knob=="stomo":
			return [ "tomo_all", "tomo_auto_only", "tomo_cross_only" ]
		elif knob=="sreduced":
			return ["none", "reduced_K10"]
		elif knob=="q_mag_size":
			print("Positive float")
			return None
		else:
			raise ValueError("{0} is not a tunable setting!".format(knob)) 




##########################################
##########Nicaea class####################
##########################################

class Nicaea(w0waCDM):

	"""
	Main class handler for the python bindings of the NICAEA cosmological code, written by M. Kilbinger & collaborators

	"""

	def __init__(self,H0=72.0,Om0=0.26,Ode0=0.74,Ob0=0.046,w0=-1.0,wa=0.0,sigma8=0.798,ns=0.960,name=None):

		if _nicaea is None:
			raise ImportError("The Nicaea bindings were not installed, check your GSL/FFTW3 installations!")

		super(Nicaea,self).__init__(H0,Om0,Ode0,w0=w0,wa=wa,name=name)
		self.sigma8=sigma8
		self.ns=ns
		self.Ob0=Ob0

	def __repr__(self):

		astropy_string = super(Nicaea,self).__repr__()
		pieces = astropy_string.split(",")
		si8_piece = u" sigma8={0}".format(self.sigma8)
		ns_piece = u" ns={0}".format(self.ns)
		Ob0_piece = u" Ob0={0}".format(self.Ob0)

		return ",".join(pieces[:3] + [si8_piece,ns_piece,Ob0_piece] + pieces[3:])

	@classmethod
	def fromCosmology(cls,cosmo):

		"""
		Builds a Nicaea instance from one of astropy.cosmology objects, from which it inherits all the cosmological parameter values

		:param cosmo: one of astropy cosmology instances
		:type cosmo: astropy FLRW

		:returns: Nicaea instance with the cosmological parameters inherited from cosmo 

		"""

		#Get the cosmological parameter values out of the cosmo object
		H0 = cosmo.H0.value
		Om0 = cosmo.Om0
		Ode0 = cosmo.Ode0

		#Dark energy
		if hasattr(cosmo,"w0"):
			w0=cosmo.w0
		else:
			w0=-1.0

		if hasattr(cosmo,"wa"):
			wa=cosmo.wa
		else:
			wa=0.0

		#Neutrinos
		Neff = cosmo.Neff
		Onu0 = cosmo.Onu0

		#Set these manually to default
		ns = 0.960
		sigma8 = 0.800
		Ob0 = 0.046

		#Instantiate
		return cls(H0=H0,Om0=Om0,Ode0=Ode0,Ob0=Ob0,w0=w0,wa=wa,sigma8=sigma8,ns=ns)


	def convergencePowerSpectrum(self,ell,z=2.0,distribution=None,settings=None):

		"""
		Computes the convergence power spectrum for the given cosmological parameters and redshift distribution using NICAEA

		:param ell: multipole moments at which to compute the power spectrum
		:type ell: array.

		:param z: redshift bins for the sources; if a single float is passed, single redshift is assumed
		:type z: float.

		:param distribution: redshift distribution of the sources (normalization not necessary)
		:type distribution: None or callable

		:param settings: NICAEA code settings
		:type settings: NicaeaSettings instance

		:returns: (array) computed power spectrum at the selected multipoles

		"""

		#If no settings provided, use the default ones
		if settings is None:
			settings=NicaeaSettings.default()

		#Parse redshift distribution from input
		if type(z)==np.float:
			
			nzbins = 1
			Nnz = np.array([2],dtype=np.int32)
			nofz = ["single"]
			par_nz = np.array([z,z])

		elif type(z)==np.ndarray:
			raise TypeError("Not implemented yet!")
		else:
			raise TypeError("Redshift format not recognized!")
		
		
		return _nicaea.shearPowerSpectrum(self.Om0,self.Ode0,self.w0,self.wa,self.H0.value/100.0,self.Ob0,self.Onu0,self.Neff,self.sigma8,self.ns,nzbins,ell,Nnz,nofz,par_nz,settings)


