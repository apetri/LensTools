import os,glob

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


#####################################################
##############SimulationModel class##################
#####################################################

class SimulationModel(object):

	"""
	Class handler of a weak lensing simulation model, defined by a set of cosmological parameters

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
				os.mkdir(d)
				print("[+] {0} created".format(d))


	def newCollection(self,box_size=240.0*u.Mpc,nside=512):

		"""
		Instantiate new simulation with the specified settings

		"""

		newSimulation = SimulationCollection(self.cosmology,self.environment,self.parameters)
		newSimulation.box_size = box_size
		newSimulation.nside = nside

		#Build the geometry_id
		newSimulation.geometry_id = "{0}b{1}".format(nside,int(box_size.to(self.Mpc_over_h).value))

		#Make the corresponding directory if not already present
		newSimulation.home_subdir = os.path.join(newSimulation.environment.home,newSimulation.cosmo_id,newSimulation.geometry_id)
		newSimulation.storage_subdir = os.path.join(newSimulation.environment.storage,newSimulation.cosmo_id,newSimulation.geometry_id)

		for d in [newSimulation.home_subdir,newSimulation.storage_subdir]:
			if not os.path.isdir(d):
				os.mkdir(d)
				print("[+] {0} created".format(d))


		return newSimulation


##########################################################
##############SimulationCollection class##################
##########################################################

class SimulationCollection(SimulationModel):


	"""
	Class handler of a collection of simulations that share model parameters

	"""

	def newCollection(self):
		raise

	def newInitialCondition(self,seed=0):
		
		#Check if there are already generated initial conditions in there
		ics_present = glob.glob(os.path.join(self.storage_subdir,"ic*"))
		new_ic_index = len(ics_present) + 1

		#Generate the new initial condition
		newIC = SimulationIC(self.cosmology,self.environment,self.parameters)

		#These are inherited from before
		newIC.box_size = self.box_size
		newIC.side = self.nside
		newIC.geometry_id = self.geometry_id

		#And these are specific to the new initial condition
		newIC.seed = seed
		newIC.storage_subdir = os.path.join(self.storage_subdir,"ic{0}".format(new_ic_index))
		os.mkdir(newIC.storage_subdir)
		print("[+] {0} created".format(newIC.storage_subdir))

		#Make new file with the number of the seed
		seedfile = open(os.path.join(newIC.storage_subdir,"seed"+str(seed)),"w")
		seedfile.close()

		return newIC



##########################################################
##############SimulationIC class##########################
##########################################################

class SimulationIC(SimulationCollection):

	"""
	Class handler of a simulation with a defined initial condition

	"""


