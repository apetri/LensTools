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
		return cls()


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


	def PowerSpectrum(self,ell,z,distribution=None,tomography=False,settings=NicaeaSettings.default()):

		nzbins=1;
		Nnz = np.array([3],dtype=np.int32)
		nofz = ["hist"]
		par_nz = np.array([2.0,2.01,10.0])
		_nicaea.shearPowerSpectrum(self.Om0,self.Ode0,self.w0,self.wa,self.H0.value/100.0,self.Ob0,self.Onu0,self.Neff,self.sigma8,self.ns,nzbins,ell,Nnz,nofz,par_nz,settings)


