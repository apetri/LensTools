import os
from astropy.cosmology import FLRW,WMAP9


from .environment import EnvironmentSettings

name2attr = dict()
name2attr["Om"] = "Om0"
name2attr["Ol"] = "Ode0"
name2attr["w"] = "w0"
name2attr["wa"] = "wa"
name2attr["h"] = "h"
name2attr["Ob"] = "Ob0"
name2attr["si"] = "sigma8"
name2attr["ns"] = "ns"


################################################
##############Simulation class##################
################################################

class Simulation(object):

	"""
	Class handler of a weak lensing simulation, defined by a set of cosmological parameters

	"""

	def __init__(self,cosmology=WMAP9,environment=None,parameters=["Om","Ol","w","ns","si"]):

		"""
		Set the base for the simulation

		:param cosmology: cosmological model to simulate
		:type cosmology: FLRW

		:param environment: environment settings of the current machine
		:type environment: EnvironmentSettings

		:param parameters: cosmological parameters to keep track of
		:type parameters: list.

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
		self.parameters = parameters

		#Build the cosmo_id
		self.cosmo_id = "_".join([ "{0}{1:.3f}".format(p,getattr(self.cosmology,name2attr[p])) for p in parameters if (hasattr(self.cosmology,name2attr[p]) and getattr(self.cosmology,name2attr[p]) is not None)])
