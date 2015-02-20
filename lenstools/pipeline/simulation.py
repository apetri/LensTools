import os

import astropy.units as u
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

		#Define the scaled unit length for convenience
		self.kpc_over_h = u.def_unit("kpc/h",u.kpc/self.cosmology.h)
		self.Mpc_over_h = u.def_unit("Mpc/h",u.Mpc/self.cosmology.h)

		#Build the cosmo_id
		self.cosmo_id = "_".join([ "{0}{1:.3f}".format(p,getattr(self.cosmology,name2attr[p])) for p in parameters if (hasattr(self.cosmology,name2attr[p]) and getattr(self.cosmology,name2attr[p]) is not None)])

		#Create directories accordingly
		self.home_subdir = os.path.join(self.environment.home,self.cosmo_id)
		self.storage_subdir = os.path.join(self.environment.storage,self.cosmo_id)

		for d in [self.home_subdir,self.storage_subdir]:
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

		#Build the geometry_id
		newSimulation.geometry_id = "{0}b{1}".format(settings.ngenic.nside,int(settings.ngenic.box_size.to(self.Mpc_over_h).value))

		#Make the corresponding directory if not already present
		newSimulation.home_subdir = os.path.join(newSimulation.environment.home,newSimulation.cosmo_id,newSimulation.geometry_id)
		newSimulation.storage_subdir = os.path.join(newSimulation.environment.storage,newSimulation.cosmo_id,newSimulation.geometry_id)

		for d in [newSimulation.home_subdir,newSimulation.storage_subdir]:
			if not os.path.isdir(d):
				print("[+] {0} created".format(d))
				os.mkdir(d)


		return newSimulation



########################################################
##############SimulationSettings class##################
########################################################

class SimulationSettings(object):

	"""
	Class handler of the simulation settings

	"""

	def __init__(self,ngenic=None,gadget=None,planes=None,rayTracing=None):

		self.ngenic = ngenic
		self.gadget = gadget
		self.planes = planes
		self.rayTracing = rayTracing


####################################################
##############NgenICSettings class##################
####################################################

class NGenICSettings(object):

	"""
	Class handler of N-GenIC settings

	"""

	def __init__(self):

		self.box_size = 240.0*u.Mpc
		self.nside = 512
		self.seed = 0


