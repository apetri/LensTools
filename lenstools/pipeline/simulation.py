import os
from astropy.cosmology import FLRW,WMAP9


from .environment import EnvironmentSettings


################################################
##############Simulation class##################
################################################

class Simulation(object):

	"""
	Class handler of a weak lensing simulation, defined by a set of cosmological parameters

	"""

	def __init__(self,cosmology=WMAP9,environment=None):

		"""
		Set the base for the simulation

		:param cosmology: cosmological model to simulate
		:type cosmology: FLRW

		:param environment: environment settings of the current machine
		:type environment: EnvironmentSettings

		"""

		#Safety checks
		assert isinstance(cosmology,FLRW)

		if environment is None:
			self.environment = EnvironmentSettings()
		else:
			assert isinstance(environment,EnvironmentSettings)
			self.environment = environment

		self.cosmology = cosmology
		self.environment = environment

		#Build the cosmo_id
		self.cosmo_id = "h{0:.3f}_Om{1:.3f}_Ol{2:.3f}".format(self.cosmology.h,self.cosmology.Om0,self.cosmology.Ode0)
		
		if (hasattr(self.cosmology,"Ob0")) and (getattr(self.cosmology,"Ob0") is not None):
			self.cosmo_id += "_Ob{0:.3f}".format(getattr(self.cosmology,"Ob0"))

		if (hasattr(self.cosmology,"w0")) and (getattr(self.cosmology,"w0") is not None):
			self.cosmo_id += "_w{0:.3f}".format(getattr(self.cosmology,"w0"))

		if (hasattr(self.cosmology,"wa")) and (getattr(self.cosmology,"wa") is not None):
			self.cosmo_id += "_wa{0:.3f}".format(getattr(self.cosmology,"wa"))

		if (hasattr(self.cosmology,"ns")) and (getattr(self.cosmology,"ns") is not None):
			self.cosmo_id += "_ns{0:.3f}".format(getattr(self.cosmology,"ns"))

		if (hasattr(self.cosmology,"sigma8")) and (getattr(self.cosmology,"sigma8") is not None):
			self.cosmo_id += "_si{0:.3f}".format(getattr(self.cosmology,"sigma8"))
