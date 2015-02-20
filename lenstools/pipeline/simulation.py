import os
from astropy.cosmology import FLRW,WMAP9

from .environment import EnvironmentSettings
from ..simulations import Gadget2Settings 

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
		self.parameters = parameters

		#Build the cosmo_id
		self.cosmo_id = "_".join([ "{0}{1:.3f}".format(p,getattr(self.cosmology,name2attr[p])) for p in parameters if (hasattr(self.cosmology,name2attr[p]) and getattr(self.cosmology,name2attr[p]) is not None)])

		#Create directories accordingly
		home_subdir = os.path.join(self.environment.home,self.cosmo_id)
		storage_subdir = os.path.join(self.environment.storage,self.cosmo_id)

		for d in [home_subdir,storage_subdir]:
			if not os.path.isdir(d):
				print("[+] {0} created".format(d))
				os.mkdir(d)


	def new(self,settings):

		"""
		Instantiate new simulation with the specified settings

		:param settings: settings of the new simulation
		:type settings: SimulationSettings

		"""

		assert isinstance(settings,SimulationSettings)
		newSimulation = self.__class__(self.cosmology,self.environment,self.parameters)
		newSimulation.settings = settings

		return newSimulation



########################################################
##############SimulationSettings class##################
########################################################

class SimulationSettings(object):

	"""
	Class handler of the simulation settings

	"""

	def __init__(self,gadget=None,planes=None,rayTracing=None):

		self.gadget = gadget
		self.planes = planes
		self.rayTracing = rayTracing


