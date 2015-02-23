import os,glob
import re

import astropy.units as u
from astropy.cosmology import FLRW,WMAP9


from .environment import EnvironmentSettings
from ..simulations import Gadget2Settings,Nicaea 

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

	@classmethod
	def available(cls,environment):

		"""
		Lists all currently available models in the home and storage directories

		:param environment: environment settings of the current machine
		:type environment: EnvironmentSettings

		:returns: list.
		:rtype: SimulationModel

		"""

		models = list()

		#Useful regular expression
		parmatch = re.compile(r"([a-zA-Z]+)([0-9.-]+)")

		if not(os.path.isdir(environment.home)) or not(os.path.isdir(environment.storage)):
			return models

		#Check all available models in the home directory
		dirnames = [ os.path.basename(n) for n in glob.glob(os.path.join(environment.home,"*")) ]

		#Decide which of the directories actually correspond to cosmological models
		for dirname in dirnames:
			
			parameters_dict = dict()
			parameters_list = list()
			parameters = dirname.split("_")
			
			for parameter in parameters:
				
				try:
					par,val = parmatch.match(parameter).groups()
				except TypeError:
					pass
				
				parameters_list.append(par)
				parameters_dict[name2attr[par]] = float(val)

			try:
				cosmoModel = Nicaea(**parameters_dict)
				models.append(cls(cosmology=cosmoModel,environment=environment,parameters=parameters_list))
			except TypeError:
				pass

		#Return the list with the available models
		return models

	def __repr__(self):

		representation_parameters = []
		for p in self.parameters:
			representation_parameters.append("{0}={1:.3f}".format(p,getattr(self.cosmology,name2attr[p])))

		return "<"+ " , ".join(representation_parameters) + ">"

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

		newSimulation = SimulationCollection(self.cosmology,self.environment,self.parameters,box_size,nside)

		#Make the corresponding directory if not already present
		for d in [newSimulation.home_subdir,newSimulation.storage_subdir]:
			if not os.path.isdir(d):
				os.mkdir(d)
				print("[+] {0} created".format(d))


		return newSimulation


	@property
	def collections(self):

		"""
		Lists all the available collections for a model

		:returns: list.
		:rtype: SimulationCollection

		"""

		collection_names = [ os.path.basename(d) for d in glob.glob(os.path.join(self.home_subdir,"*")) ]
		collection_list = list()

		for name in collection_names:
			
			try:
				nside,box = name.split("b")
				collection_list.append(SimulationCollection(self.cosmology,self.environment,self.parameters,box_size=float(box)*self.Mpc_over_h,nside=int(nside)))
			except ValueError:
				pass

		return collection_list


##########################################################
##############SimulationCollection class##################
##########################################################

class SimulationCollection(SimulationModel):


	"""
	Class handler of a collection of simulations that share model parameters

	"""

	def __init__(self,cosmology,environment,parameters,box_size,nside):

		super(SimulationCollection,self).__init__(cosmology,environment,parameters)
		self.box_size = box_size
		self.nside = nside

		#Build the geometry_id
		self.geometry_id = "{0}b{1}".format(nside,int(box_size.to(self.Mpc_over_h).value))

		#Build the directory names
		self.home_subdir = os.path.join(self.environment.home,self.cosmo_id,self.geometry_id)
		self.storage_subdir = os.path.join(self.environment.storage,self.cosmo_id,self.geometry_id)


	def __repr__(self):

		return super(SimulationCollection,self).__repr__() + " ; box={0},nside={1}".format(self.box_size,self.nside)


	def newCollection(self):
		raise TypeError("This method should be called on SimulationModel instances!")

	def newInitialCondition(self,seed=0):
		
		#Check if there are already generated initial conditions in there
		ics_present = glob.glob(os.path.join(self.storage_subdir,"ic*"))
		new_ic_index = len(ics_present) + 1

		#Generate the new initial condition
		newIC = SimulationIC(self.cosmology,self.environment,self.parameters,self.box_size,self.nside,new_ic_index,seed)

		#Make dedicated directory for new initial condition
		os.mkdir(newIC.storage_subdir)
		print("[+] {0} created".format(newIC.storage_subdir))

		#Make new file with the number of the seed
		seedfile = open(os.path.join(newIC.storage_subdir,"seed"+str(seed)),"w")
		seedfile.close()

		return newIC


	@property
	def ics(self):

		"""
		List the available initial conditions (or independent simulations) for the current collection

		:returns: list.
		:rtype: SimulationIC

		"""

		ic_list = list()
		ic_numbers = [ os.path.basename(d).strip("ic") for d in glob.glob(os.path.join(self.storage_subdir,"ic*")) ]

		for ic in ic_numbers:

			seed = int(os.path.basename(glob.glob(os.path.join(self.storage_subdir,"ic"+ic,"seed*"))[0]).strip("seed"))
			ic_list.append(SimulationIC(self.cosmology,self.environment,self.parameters,self.box_size,self.nside,int(ic),seed))

		return ic_list



##########################################################
##############SimulationIC class##########################
##########################################################

class SimulationIC(SimulationCollection):

	"""
	Class handler of a simulation with a defined initial condition

	"""

	def __init__(self,cosmology,environment,parameters,box_size,nside,ic_index,seed):

		super(SimulationIC,self).__init__(cosmology,environment,parameters,box_size,nside)

		#Save random seed information as attribute
		self.ic_index = ic_index
		self.seed = seed

		#Save storage sub-directory name
		self.storage_subdir = os.path.join(self.storage_subdir,"ic{0}".format(ic_index))

	def __repr__(self):

		return super(SimulationIC,self).__repr__() + " ; ic={0},seed={1}".format(self.ic_index,self.seed)

	def newInitialCondition(self,seed):
		raise TypeError("This method should be called on SimulationCollection instances!")


