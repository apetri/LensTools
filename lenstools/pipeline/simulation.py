from __future__ import division,with_statement
from abc import ABCMeta,abstractproperty,abstractmethod

import sys,os
import re
import tarfile
import json
import itertools

if sys.version_info.major>=3:
	from io import StringIO
else:
	from StringIO import StringIO

import numpy as np
import astropy.units as u
from astropy.cosmology import z_at_value

from .. import configuration
from ..utils.configuration import LensToolsCosmology

from .remote import SystemHandler,LocalGit
from .settings import *

from .deploy import JobHandler

from ..simulations.camb import CAMBTransferFromPower
from ..simulations import Gadget2SnapshotDE
from ..simulations.raytracing import PotentialPlane

#####################################################
############Parse cosmology from string##############
#####################################################

def string2cosmo(s,name2attr):

	parmatch = re.compile(r"([a-zA-Z]+)([0-9.-]+)")

	parameters_dict = dict()
	parameters_list = list()
	parameters = s.split("_")

	for parameter in parameters:
				
		try:
			par,val = parmatch.match(parameter).groups()
		except AttributeError:
			return None
				
		parameters_list.append(par)

		try:
			
			#The Hubble parameter needs particular attention
			if par=="h":
				parameters_dict["H0"] = 100.0*float(val)
			else: 
				parameters_dict[name2attr[par]] = float(val)
		
		except (ValueError,KeyError):
			return None

	try:
		cosmoModel = LensToolsCosmology(**parameters_dict)
	except TypeError:
		return None

	return cosmoModel,parameters_list


###################################################################
############CAMB --> NGenIC power spectrum conversion##############
###################################################################

def _camb2ngenic(k,P):

	lgk = np.log(k) / np.log(10)
	lgP = np.log((k**3)*P/(2.0*(np.pi**2))) / np.log(10)

	return lgk,lgP

##############################################
##############InfoDict class##################
##############################################

class InfoDict(object):

	def __init__(self,batch,**kwargs):

		if isinstance(batch,SimulationBatch):
			self.batch = batch
		elif isinstance(batch,EnvironmentSettings):
			self.batch = SimulationBatch(batch,**kwargs)
		else:
			raise TypeError("batch type not recognized!")

		self.dictionary = self.batch.info

	def __enter__(self):
		return self

	def __exit__(self,type,value,tb):
		pass

	def update(self):
		dictionary_file = os.path.join(self.batch.home_subdir,self.batch.environment.json_tree_file)
		with self.batch.syshandler.open(dictionary_file,"w") as fp:
			fp.write(json.dumps(self.dictionary))

#############################################################
##############Convenient resource retrieval##################
#############################################################

def _retrieve(self,s):

	#Break down search string into single components
	search = "mcrpMCS"
	nsteps = len(filter(lambda c:c in search,s))

	if not nsteps:
		return None

	s_pieces = re.match(r"([{0}][0-9]+)".format(search)*nsteps,s).groups()

	#Walk to the pieces and return the corresponding resource
	current = self
	for p in s_pieces:

		#Split into resoure and index number
		resource,index = re.match(r"([{0}])([0-9]+)".format(search),p).groups()

		#model
		if resource=="m":
			current = current.models[int(index)]

		#collection
		if resource=="c":
			current = current.collections[int(index)]

		#realization
		if resource=="r":
			current = current.realizations[int(index)]

		#plane set
		if resource=="p":
			current = current.planesets[int(index)]

		#map set
		if resource=="M":
			current = current.mapsets[int(index)]

		#catalog
		if resource=="C":
			current = current.catalogs[int(index)]

		#sub-catalog
		if resource=="S":
			current = current.subcatalogs[int(index)]

	#Return resource to user
	return current

#####################################################
##############SimulationBatch class##################
#####################################################

class SimulationBatch(object):

	"""
	Class handler of a batch of weak lensing simulations that share the same environment settings

	"""

	#Keep track of which batches are already loaded in memory (give the class a singleton--like aspect)
	_in_memory = dict()

	@classmethod
	def current(cls,name="environment.ini",syshandler=configuration.syshandler,indicize=False):

		"""
		This method looks in the current directory and looks for a configuration file named "environment.ini"; if it finds one, it returns a SimulationBatch instance that corresponds to the one pointed to by "environment.ini" (default)

		:param name: name of the INI file with the environment settings, defaults to 'environment.ini'
		:type name: str.

		:param syshandler: system handler that allows to override the methods used to create directories and do I/O from files, must implement the abstract type SystemHandler
		:type syshandler: SystemHandler

		:returns: Simulation batch pointed to by "environment.ini", or None
		:rtype: SimulationBatch

		"""

		if not(os.path.exists(name)):
			return None

		env = EnvironmentSettings.read(name)
		return cls(env,syshandler,indicize)

	@property 
	def home_subdir(self):
		return self.environment.home

	@property
	def storage_subdir(self):
		return self.environment.storage

	@property
	def home(self):
		return self.home_subdir

	@property
	def storage(self):
		return self.storage_subdir
		

	def __init__(self,environment,syshandler=configuration.syshandler,indicize=False):

		"""
		Gets the handler instance of a batch of simulations residing in the provided environment

		:param environment: environment settings
		:type environment: EnvironmentSettings

		:param syshandler: system handler that allows to override the methods used to create directories and do I/O from files, must implement the abstract type SystemHandler
		:type syshandler: SystemHandler

		"""

		#Type check
		assert isinstance(environment,EnvironmentSettings)
		assert isinstance(syshandler,SystemHandler)

		if environment.home in self.__class__._in_memory.keys():
			return
		
		self.environment = environment
		self.syshandler = syshandler

		#Create directories if they do not exist yet
		if not self.syshandler.isbatch(environment.home):
			
			self.syshandler.init(environment.home)
			print("[+] {0} created on {1}".format(environment.home,self.syshandler.name))

			#Create also an "environment.ini" file that provides easy access to the current simulation batch from Home
			with self.syshandler.open(os.path.join(environment.home,"environment.ini"),"w") as envfile:
				envfile.write("[EnvironmentSettings]\n\n")
				envfile.write("home = {0}\n".format(os.path.abspath(environment.home)))
				envfile.write("storage = {0}\n\n".format(os.path.abspath(environment.storage)))


		if not self.syshandler.exists(environment.storage):
			self.syshandler.mkdir(environment.storage)
			print("[+] {0} created on {1}".format(environment.storage,self.syshandler.name))

		#Indicize the simulation products
		if indicize:
			self.update_changes()

		#Keep track of this simulation batch
		self.__class__._in_memory[self.home_subdir] = self


	def __new__(cls,environment,*args,**kwargs):

		if environment.home in cls._in_memory:
			return cls._in_memory[environment.home]
		else:
			return object.__new__(cls)

	

	def update_changes(self):
		with InfoDict(self) as info:
			info.update()

	##############################################################################################################################

	#Convenient resource retrieval
	def __getitem__(self,s):
		return _retrieve(self,s)

	##############################################################################################################################

	def commit(self,message):

		"""
		If the Simulation Batch is put under version control in a git repository, this method commits the newly added models,collections,realizations or map/plane sets

		:param message: commit message
		:type message: str.

		"""

		if isinstance(self.syshandler,LocalGit):
			self.syshandler.repository.index.commit(message)
		else:
			raise TypeError("The system handler must be an instance of LocalGit to use this method!")

	##############################################################################################################################

	@property
	def available(self):

		"""
		Lists all currently available models in the home and storage directories

		:returns: list.
		:rtype: SimulationModel

		"""

		models = list()

		#Check all available models in the home directory
		if self.syshandler.exists(self.infofile):
			dirnames = self.info.keys()
		else:
			dirnames = [ os.path.basename(n) for n in self.syshandler.glob(os.path.join(self.environment.home,"*")) ]

		#Cycle over directory names
		for d in dirnames:
			model = self.getModel(d)
			if model is not None:
				models.append(model)

		return models

	#Alias for available models
	@property
	def models(self):
		return self.available

	##############################################################################################################################
	@property
	def infofile(self):
		return os.path.join(self.home_subdir,self.environment.json_tree_file)

	@property
	def info(self):

		"""
		Returns summary info of the simulation batch corresponding to the current environment

		:returns: info in dictionary format

		"""

		#See if the info is already available
		if hasattr(self,"_info"):
			return self._info

		#Load info from tree file if available
		tree_file_path = os.path.join(self.home_subdir,self.environment.json_tree_file)
		if self.syshandler.exists(tree_file_path):
			with self.syshandler.open(tree_file_path,"r") as fp:
				self._info = json.loads(fp.read())
				return self._info

		#Information will be returned in dictionary format
		info_dict = dict()

		#Start with the available models
		available_models = self.available
		for model in available_models:
			info_dict[model.cosmo_id] = dict()

			#Follow with the collections 
			for collection in model.collections:
				info_dict[model.cosmo_id][collection.geometry_id] = dict()
				info_dict[model.cosmo_id][collection.geometry_id]["nbody"] = dict()
				info_dict[model.cosmo_id][collection.geometry_id]["map_sets"] = dict()
				info_dict[model.cosmo_id][collection.geometry_id]["catalogs"] = dict()

				#Check if there are any map sets or catalogs present
				try:
					
					with self.syshandler.open(os.path.join(collection.home_subdir,"sets.txt"),"r") as setsfile:
						
						for line in setsfile.readlines():
							if line=="":
								continue
							map_set = collection.getMapSet(line.strip("\n"))
							info_dict[model.cosmo_id][collection.geometry_id]["map_sets"][map_set.settings.directory_name] = dict()
							info_dict[model.cosmo_id][collection.geometry_id]["map_sets"][map_set.settings.directory_name]["settings"] = map_set.settings.to_dict()

				except IOError:
					pass

				
				try:
					
					with self.syshandler.open(os.path.join(collection.home_subdir,"catalogs.txt"),"r") as setsfile:
						
						for line in setsfile.readlines():
							if line=="":
								continue
							
							catalog = collection.getCatalog(line.strip("\n"))
							info_dict[model.cosmo_id][collection.geometry_id]["catalogs"][catalog.settings.directory_name] = dict() 
							info_dict[model.cosmo_id][collection.geometry_id]["catalogs"][catalog.settings.directory_name]["settings"] = catalog.settings.to_dict()

				except IOError:
					pass

				#Follow with the realizations
				for r in collection.realizations:
					
					info_dict[model.cosmo_id][collection.geometry_id]["nbody"][r.ic_index] = dict()

					try:
						info_dict[model.cosmo_id][collection.geometry_id]["nbody"][r.ic_index]["ngenic_settings"] = r.ngenic_settings.to_dict()
						info_dict[model.cosmo_id][collection.geometry_id]["nbody"][r.ic_index]["gadget_settings"] = r.gadget_settings.to_dict()
					except AttributeError:
						pass
					
					info_dict[model.cosmo_id][collection.geometry_id]["nbody"][r.ic_index]["plane_sets"] = dict()

					#Check if there are any plane sets present
					try:

						with self.syshandler.open(os.path.join(r.home_subdir,"sets.txt"),"r") as setsfile:

							for line in setsfile.readlines():
								if line=="":
									continue
								
								plane_set = r.getPlaneSet(line.strip("\n"))
								info_dict[model.cosmo_id][collection.geometry_id]["nbody"][r.ic_index]["plane_sets"][plane_set.settings.directory_name] = dict()
								info_dict[model.cosmo_id][collection.geometry_id]["nbody"][r.ic_index]["plane_sets"][plane_set.settings.directory_name]["settings"] = plane_set.settings.to_dict()

					except IOError:
						pass


		#Return to user
		return info_dict


	##############################################################################################################################################

	def list(self,resource=None,which=None,chunk_size=10,**kwargs):

		"""
		Lists the available resources in the simulation batch (collections,mapsets,etc...)

		:param resource: custom function to call on each batch.models element, must return a string. If None the list of Storage model directories is returned
		:type resource: None or callable

		:param which: extremes of the model numbers to get (if None all models are processed); if callable, filter(which,self.models) gives the models to archive
		:type which: tuple. or callable

		:param chunk_size: size of output chunk
		:type chunk_size: int.

		:param kwargs: the keyword arguments are passed to resource
		:type kwargs: dict.

		:returns: requested resources
		:rtype: list.

		"""

		#Available models
		if which is None:
			models = self.models
		elif isinstance(which,tuple):
			models = self.models.__getslice__(*which)
		else:
			models = filter(which,self.models)


		#Return chunks
		chunks = list()
		local_chunk = list()

		while True:

			#Get the model at the front
			try:
				model = models.pop(0)
			except IndexError:
				if len(local_chunk):
					chunks.append("\n".join(local_chunk))
				break

			#Extract the resource
			if resource is not None:
				local_chunk.append(resource(model,**kwargs))
			else:
				local_chunk.append(model.storage_subdir)

			#If we reached the chunk size dump and reset
			if len(local_chunk)==chunk_size:
				chunks.append("\n".join(local_chunk))
				local_chunk = list()

		#Return to user
		return chunks


	#################
	####Compress#####
	#################

	@staticmethod
	def _archive(name,chunk,mode):

		print("[+] Compressing {0} into {1}".format("-".join(chunk.split("\n")),name))

		with tarfile.open(name,mode) as tar:
			for f in chunk.split("\n"):
				tar.add(f)

	##############################################################################################################################################

	def archive(self,name,pool=None,chunk_size=1,**kwargs):

		"""
		Archives a batch available resource to a tar gzipped archive; the resource file/directory names are retrieved with the list method. The archives will be written to the simulation batch storage directory 

		:param name: name of the archive
		:type name: str.

		:param pool: MPI Pool used to parallelize compression
		:type pool: MPIPool

		:param kwargs: the keyword arguments are passed to the list method
		:type kwargs: dict.

		"""

		#Retrieve resource chunks
		resource_chunks = self.list(chunk_size=chunk_size,**kwargs)

		#Get archive names
		if type(name)==str:

			name_pieces = name.split(".")
			if name_pieces[-1]=="gz":
				mode = "w"
			else:
				mode = "w:gz"
			
			archive_names = list()

			for n,chunk in enumerate(resource_chunks):

				#Build archive name
				archive_name = name.replace(".tar","{0}.tar".format(n+1))
				archive_names.append(os.path.join(self.environment.storage,archive_name))

		elif type(name)==list:

			mode = "w:gz"
			archive_names = [ os.path.join(self.environment.storage,n) for n in name ]

		else:
			raise TypeError("name should be a string or list!")

		#Safety assert
		assert len(archive_names)==len(resource_chunks),"You should provide an archive file name for each resouce chunk!"
		
		#Call the _archive method to make the compression
		if pool is None:

			for n,chunk in enumerate(resource_chunks):
				self.__class__._archive(archive_names[n],chunk,mode)
		
		else:
			assert len(resource_chunks)==pool.size+1,"There should be one MPI task (you have {0}) for each chunk (you have {1})!".format(pool.size+1,len(resource_chunks))
			self.__class__._archive(archive_names[pool.rank],resource_chunks[pool.rank],mode)

	#################
	####Unpack#######
	#################

	@staticmethod
	def _unpack(name,path):
		print("[+] Unpacking {0} into {1}...".format(name,path))
		with tarfile.open(name,"r:gz") as tar:
			tar.extractall(path=path)

	def unpack(self,where,which=None,pool=None):

		"""
		Unpacks the compressed simulation batch products into the storage directory: the resources of each model must be contained in a file called <cosmo_id>.tar.gz

		:param where: path of the compressed resources
		:type where: str.

		:param which: extremes of the model numbers to unpack (if None all models are unpacked)
		:type which: tuple.

		:param pool: MPI Pool used to parallelize de-compression
		:type pool: MPIPool

		"""

		#Get the models to unpack 
		if which is None:
			models = self.available
		else:
			models = self.available.__getslice__(*which)


		#Unpack each model (spread computations over MPIPool if provided)
		if pool is not None:
			assert len(models)==pool.size+1,"The number of MPI processes must be equal to the number of models to unpack!"

		if pool is None:
			for model in models:
				archive_path = os.path.join(where,"{0}.tar.gz".format(model.cosmo_id))
				self._unpack(archive_path,self.environment.storage)
		else:
			archive_path = os.path.join(where,"{0}.tar.gz".format(models[pool.rank].cosmo_id))
			self._unpack(archive_path,self.environment.storage)



	##############################################################################################################################################

	####################
	####Duplicate#######
	####################

	def copyTree(self,path,syshandler=configuration.syshandler):

		"""
		Copies the current batch directory tree into a separate path

		:param path: path into which to copy the current batch directory tree
		:type path: str.

		:param syshandler: system handler (can be a remote)
		:type syshandler: SystemHandler

		"""

		#Instantiate new SimulationBatch object (home and storage will be the same)
		environment = EnvironmentSettings(home=path,storage=path)
		for key in ["cosmo_id_digits","name2attr","json_tree_file"]:
			setattr(environment,key,getattr(self.environment,key))
			
		batchCopy = SimulationBatch(environment,syshandler)

		#Walk down the directory tree and create the copied directories on the go (only if non existent already)

		#Model
		for model in self.available:

			modelCopy = batchCopy.getModel(model.cosmo_id)

			if modelCopy is None:
				modelCopy = batchCopy.newModel(model.cosmology,model.parameters)

			#Collection
			for coll in model.collections:

				collCopy = modelCopy.getCollection(box_size=coll.box_size,nside=coll.nside)

				if collCopy is None: 
					collCopy = modelCopy.newCollection(box_size=coll.box_size,nside=coll.nside)

				#Maps
				for map_set in coll.mapsets:

					map_setCopy = collCopy.getMapSet(map_set.name)

					if map_setCopy is None:
						map_setCopy = collCopy.newMapSet(map_set.settings)

				#Catalogs
				for catalog in coll.catalogs:

					catalogCopy = collCopy.getCatalog(catalog.name)

					if catalogCopy is None:
						catalogCopy = collCopy.newCatalog(catalog.settings)

				#Realizations
				for r in coll.realizations:

					rCopy = collCopy.getRealization(r.ic_index)

					if rCopy is None:
						rCopy = collCopy.newRealization(seed=r.seed)

					#Planes
					for plane_set in r.planesets:

						plane_setCopy = rCopy.getPlaneSet(plane_set.name)

						if plane_setCopy is None:
							plane_setCopy = rCopy.newPlaneSet(plane_set.settings)


		#Return handle on the copied directory tree
		return batchCopy


	##############################################################################################################################################

	def newModel(self,cosmology,parameters):

		"""
		Create a new simulation model, given a set of cosmological parameters

		:param cosmology: cosmological model to simulate
		:type cosmology: LensToolsCosmology

		:param parameters: cosmological parameters to keep track of
		:type parameters: list.

		:rtype: SimulationModel

		"""

		#cosmology needs to be of type LensToolsCosmology
		assert isinstance(cosmology,LensToolsCosmology)

		newModel = SimulationModel(cosmology=cosmology,environment=self.environment,parameters=parameters,syshandler=self.syshandler)

		for d in [newModel.home_subdir,newModel.storage_subdir]:
			if not self.syshandler.exists(d):

				self.syshandler.mkdir(d)
				print("[+] {0} created on {1}".format(d,self.syshandler.name))

				#Update dictionary if the batch is indicized
				if self.syshandler.exists(self.infofile):
					self.info[newModel.cosmo_id] = dict()

			else:
				print("[-] Model {0} already exists!".format(newModel.cosmo_id))		

		#Return to user
		return newModel

	###########################################################################################################################################

	def getModel(self,cosmo_id):

		"""
		Instantiate a SimulationModel object corresponding to the cosmo_id provided

		:param cosmo_id: cosmo_id of the model
		:type cosmo_id: str.

		:rtype: SimulationModel

		"""
		
		if not(self.syshandler.exists(os.path.join(self.environment.home,cosmo_id))) or not(self.syshandler.exists(os.path.join(self.environment.home,cosmo_id))):
			return None

		#Parse the cosmological model from the directory name
		cosmo_parsed = string2cosmo(cosmo_id,self.environment.name2attr)

		#Return the SimulationModel instance
		if cosmo_parsed is not None:
			return SimulationModel(cosmology=cosmo_parsed[0],environment=self.environment,parameters=cosmo_parsed[1],syshandler=self.syshandler)
		else:
			return None

	###########################################################################################################################################
	##########################################Job submission scripts###########################################################################
	###########################################################################################################################################

	def writeCAMBSubmission(self,realization_list,job_settings,job_handler,config_file="camb.param",chunks=1,**kwargs):

		"""
		Writes CAMB submission script

		:param realization_list: list of ics to generate in the form "cosmo_id|geometry_id"
		:type realization_list: list. of str.

		:param chunks: number of independent jobs in which to split the submission (one script per job will be written)
		:type chunks: int. 

		:param job_settings: settings for the job (resources, etc...)
		:type job_settings: JobSettings

		:param job_handler: handler of the cluster specific features (job scheduler, architecture, etc...)
		:type job_handler: JobHandler

		:param config_file: name of the CAMB configuration file
		:type config_file: str.

		:param kwargs: you can set one_script=True to include all the executables sequentially in a single script
		:type kwargs: dict.

		"""

		#Type safety check
		assert isinstance(job_settings,JobSettings)
		assert isinstance(job_handler,JobHandler)
		assert len(realization_list)%chunks==0,"Perfect load balancing enforced, each job will process the same number of realizations!"

		#Check if we need to collapse everyting in one script
		if "one_script" in kwargs.keys():
			one_script = kwargs["one_script"]
		else:
			one_script = False

		#It's better to run CAMB from the directory where the executable resides
		job_handler.cluster_specs.execution_preamble = "cd {0}".format(os.path.dirname(job_settings.path_to_executable))

		#This limit is enforced
		assert job_settings.cores_per_simulation*chunks==len(realization_list),"cores_per_simulation x chunks should be equal to the total number of models!"
		job_settings.num_cores = job_settings.cores_per_simulation

		#Create the dedicated Job and Logs directories if not existent already
		for d in [os.path.join(self.environment.home,"Jobs"),os.path.join(self.environment.home,"Logs")]:
			
			if not(self.syshandler.exists(d)):
				self.syshandler.mkdir(d)
				print("[+] {0} created on {1}".format(d,self.syshandler.name))

		#Split realizations between independent jobs
		realizations_per_chunk = len(realization_list)//chunks

		for c in range(chunks):
		
			#Arguments for the executable
			exec_args = list()

			for realization in realization_list[realizations_per_chunk*c:realizations_per_chunk*(c+1)]:
			
				#Separate the cosmo_id,geometry_id,realization number
				cosmo_id,geometry_id = realization.split("|")

				#Get the corresponding SimulationXXX instances
				model = self.getModel(cosmo_id)

				nside,box_size = geometry_id.split("b")
				nside = int(nside)
				box_size = float(box_size) * model.Mpc_over_h
				collection = model.getCollection(box_size=box_size,nside=nside)

				parameter_file = os.path.join(collection.home_subdir,config_file)
				if not(self.syshandler.exists(parameter_file)):
					raise IOError("CAMB parameter file at {0} does not exist yet!".format(parameter_file))

				exec_args.append(parameter_file)

			executable = job_settings.path_to_executable + " " + " ".join(exec_args)

			#Write the script
			script_filename = os.path.join(self.environment.home,"Jobs",job_settings.job_script_file)
			if not one_script:
				script_filename_split = script_filename.split(".")
				script_filename_split[-2] += "{0}".format(c+1)
				script_filename = ".".join(script_filename_split)

			#Override settings to make stdout and stderr go in the right places
			job_settings.redirect_stdout = os.path.join(self.environment.home,"Logs",job_settings.redirect_stdout)
			job_settings.redirect_stderr = os.path.join(self.environment.home,"Logs",job_settings.redirect_stderr)

			if (not one_script) or (not c):
			
				#Inform user where logs will be directed
				print("[+] Stdout will be directed to {0}".format(job_settings.redirect_stdout))
				print("[+] Stderr will be directed to {0}".format(job_settings.redirect_stderr))

				with self.syshandler.open(script_filename,"w") as scriptfile:
					scriptfile.write(job_handler.writePreamble(job_settings))
			
			with self.syshandler.open(script_filename,"a") as scriptfile:
				scriptfile.write(job_handler.writeExecution([executable],[job_settings.num_cores],job_settings))

			#Log to user and return
			if (not one_script) or (not c):	
				print("[+] {0} written on {1}".format(script_filename,self.syshandler.name))


	############################################################################################################################################

	def writeNGenICSubmission(self,realization_list,job_settings,job_handler,config_file="ngenic.param",chunks=1,**kwargs):

		"""
		Writes NGenIC submission script

		:param realization_list: list of ics to generate in the form "cosmo_id|geometry_id|icN"
		:type realization_list: list. of str.

		:param chunks: number of independent jobs in which to split the submission (one script per job will be written)
		:type chunks: int. 

		:param job_settings: settings for the job (resources, etc...)
		:type job_settings: JobSettings

		:param job_handler: handler of the cluster specific features (job scheduler, architecture, etc...)
		:type job_handler: JobHandler

		:param config_file: name of the NGenIC config file
		:type config_file: str.

		:param kwargs: you can set one_script=True to include all the executables sequentially in a single script
		:type kwargs: dict.

		"""

		#Type safety check
		assert isinstance(job_settings,JobSettings)
		assert isinstance(job_handler,JobHandler)
		assert len(realization_list)%chunks==0,"Perfect load balancing enforced, each job will process the same number of realizations!"

		#Check if we need to collapse everyting in one script
		if "one_script" in kwargs.keys():
			one_script = kwargs["one_script"]
		else:
			one_script = False

		#Create the dedicated Job and Logs directories if not existent already
		for d in [os.path.join(self.environment.home,"Jobs"),os.path.join(self.environment.home,"Logs")]:
			
			if not(self.syshandler.exists(d)):
				self.syshandler.mkdir(d)
				print("[+] {0} created on {1}".format(d,self.syshandler.name))

		#Split realizations between independent jobs
		realizations_per_chunk = len(realization_list)//chunks

		for c in range(chunks):
		
			#Arguments for executable
			exec_args = list()

			for realization in realization_list[realizations_per_chunk*c:realizations_per_chunk*(c+1)]:
			
				#Separate the cosmo_id,geometry_id,realization number
				cosmo_id,geometry_id,ic_number = realization.split("|")

				#Get the corresponding SimulationXXX instances
				model = self.getModel(cosmo_id)

				nside,box_size = geometry_id.split("b")
				nside = int(nside)
				box_size = float(box_size) * model.Mpc_over_h
				collection = model.getCollection(box_size=box_size,nside=nside)

				r = collection.getRealization(int(ic_number.strip("ic")))

				parameter_file = os.path.join(r.home_subdir,config_file)
				if not(self.syshandler.exists(parameter_file)):
					raise IOError("NGenIC parameter file at {0} does not exist yet!".format(parameter_file))

				exec_args.append(parameter_file)
			
			#Executable
			executables = [ job_settings.path_to_executable + " " + " ".join(exec_args) ]

			#Write the script
			script_filename = os.path.join(self.environment.home,"Jobs",job_settings.job_script_file)
			if not one_script:
				script_filename_split = script_filename.split(".")
				script_filename_split[-2] += "{0}".format(c+1)
				script_filename = ".".join(script_filename_split)

			#Override settings to make stdout and stderr go in the right places
			job_settings.redirect_stdout = os.path.join(self.environment.home,"Logs",job_settings.redirect_stdout)
			job_settings.redirect_stderr = os.path.join(self.environment.home,"Logs",job_settings.redirect_stderr)

			#Override settings
			job_settings.num_cores = job_settings.cores_per_simulation

			if (not one_script) or (not c):
			
				#Inform user where logs will be directed
				print("[+] Stdout will be directed to {0}".format(job_settings.redirect_stdout))
				print("[+] Stderr will be directed to {0}".format(job_settings.redirect_stderr))

				with self.syshandler.open(script_filename,"w") as scriptfile:
					scriptfile.write(job_handler.writePreamble(job_settings))
			
			with self.syshandler.open(script_filename,"a") as scriptfile:
				scriptfile.write(job_handler.writeExecution(executables,[job_settings.num_cores]*len(executables),job_settings))

			#Log to user and return
			if (not one_script) or (not c):	
				print("[+] {0} written on {1}".format(script_filename,self.syshandler.name))


	############################################################################################################################################


	def writeNbodySubmission(self,realization_list,job_settings,job_handler,config_file="gadget2.param",chunks=1,**kwargs):
		
		"""
		Writes N-body simulation submission script

		:param realization_list: list of ics to generate in the form "cosmo_id|geometry_id|icN"
		:type realization_list: list. of str.

		:param chunks: number of independent jobs in which to split the submission (one script per job will be written)
		:type chunks: int. 

		:param job_settings: settings for the job (resources, etc...)
		:type job_settings: JobSettings

		:param job_handler: handler of the cluster specific features (job scheduler, architecture, etc...)
		:type job_handler: JobHandler

		:param config_file: name of the Nbody configuration file
		:type config_file: str.

		:param kwargs: you can set one_script=True to include all the executables sequentially in a single script
		:type kwargs: dict.

		"""

		#Type safety check
		assert isinstance(job_settings,JobSettings)
		assert isinstance(job_handler,JobHandler)
		assert len(realization_list)%chunks==0,"Perfect load balancing enforced, each job will process the same number of realizations!"

		#Check if we need to collapse everyting in one script
		if "one_script" in kwargs.keys():
			one_script = kwargs["one_script"]
		else:
			one_script = False

		#Create the dedicated Job and Logs directories if not existent already
		for d in [os.path.join(self.environment.home,"Jobs"),os.path.join(self.environment.home,"Logs")]:
			
			if not(self.syshandler.exists(d)):
				self.syshandler.mkdir(d)
				print("[+] {0} created on {1}".format(d,self.syshandler.name))

		#Split realizations between independent jobs
		realizations_per_chunk = len(realization_list)//chunks

		#Override job settings to make sure requested resources are enough
		job_settings.num_cores = job_settings.cores_per_simulation*realizations_per_chunk

		for c in range(chunks):

			#Find the realizations this chunk must process
			realizations_in_chunk = realization_list[realizations_per_chunk*c:realizations_per_chunk*(c+1)]

			#Find in how many executables we should split the realizations in this chunk
			num_executables = len(realizations_in_chunk)

			#Find how many realizations per executable
			executables = list()
			cores = [job_settings.cores_per_simulation] * num_executables
			
			for e in range(num_executables):
			
				#Separate the cosmo_id,geometry_id,realization number
				cosmo_id,geometry_id,ic_number = realizations_in_chunk[e].split("|")

				#Get the corresponding SimulationXXX instances
				model = self.getModel(cosmo_id)

				nside,box_size = geometry_id.split("b")
				nside = int(nside)
				box_size = float(box_size) * model.Mpc_over_h
				collection = model.getCollection(box_size=box_size,nside=nside)

				r = collection.getRealization(int(ic_number.strip("ic")))

				parameter_file = os.path.join(r.home_subdir,config_file)
				if not(self.syshandler.exists(parameter_file)):
					raise IOError("Gadget2 parameter file at {0} does not exist yet!".format(parameter_file))

				if issubclass(configuration.snapshot_handler,Gadget2SnapshotDE):
					executables.append(job_settings.path_to_executable + " " + "{0} {1} {2}".format(1,job_settings.cores_per_simulation,parameter_file))
				else:
					executables.append(job_settings.path_to_executable + " " + "{0}".format(parameter_file))

			#Write the script
			script_filename = os.path.join(self.environment.home,"Jobs",job_settings.job_script_file)

			if not one_script:
				script_filename_split = script_filename.split(".")
				script_filename_split[-2] += "{0}".format(c+1)
				script_filename = ".".join(script_filename_split)

			#Override settings to make stdout and stderr go in the right places
			job_settings.redirect_stdout = os.path.join(self.environment.home,"Logs",job_settings.redirect_stdout)
			job_settings.redirect_stderr = os.path.join(self.environment.home,"Logs",job_settings.redirect_stderr)

			if (not one_script) or (not c):
			
				#Inform user where logs will be directed
				print("[+] Stdout will be directed to {0}".format(job_settings.redirect_stdout))
				print("[+] Stderr will be directed to {0}".format(job_settings.redirect_stderr))

				with self.syshandler.open(script_filename,"w") as scriptfile:
					scriptfile.write(job_handler.writePreamble(job_settings))
			
			with self.syshandler.open(script_filename,"a") as scriptfile:
				scriptfile.write(job_handler.writeExecution(executables,cores,job_settings))

			#Log to user and return
			if (not one_script) or (not c):	
				print("[+] {0} written on {1}".format(script_filename,self.syshandler.name))

	############################################################################################################################################

	def writePlaneSubmission(self,realization_list,job_settings,job_handler,chunks=1,**kwargs):
		
		"""
		Writes lens plane generation submission script

		:param realization_list: list of ics to generate in the form "cosmo_id|geometry_id|icN"
		:type realization_list: list. of str.

		:param chunks: number of independent jobs in which to split the submission (one script per job will be written)
		:type chunks: int. 

		:param job_settings: settings for the job (resources, etc...)
		:type job_settings: JobSettings

		:param job_handler: handler of the cluster specific features (job scheduler, architecture, etc...)
		:type job_handler: JobHandler

		:param kwargs: keyword arguments accepted are "environment_file" to specify the environment settings for the current batch, "plane_config_file" to specify the lensing option for plane generation script. Additionally you can set one_script=True to include all the executables sequentially in a single script
		:type kwargs: dict.

		"""

		#Type safety check
		assert isinstance(job_settings,JobSettings)
		assert isinstance(job_handler,JobHandler)
		assert len(realization_list)%chunks==0,"Perfect load balancing enforced, each job will process the same number of realizations!"

		#Check if we need to collapse everyting in one script
		if "one_script" in kwargs.keys():
			one_script = kwargs["one_script"]
		else:
			one_script = False

		#Create the dedicated Job and Logs directories if not existent already
		for d in [os.path.join(self.environment.home,"Jobs"),os.path.join(self.environment.home,"Logs")]:
			
			if not(self.syshandler.exists(d)):
				self.syshandler.mkdir(d)
				print("[+] {0} created on {1}".format(d,self.syshandler.name))

		#Split realizations between independent jobs
		realizations_per_chunk = len(realization_list)//chunks

		#Override job settings to make sure requested resources are enough
		job_settings.num_cores = job_settings.cores_per_simulation*realizations_per_chunk

		for c in range(chunks):

			#Find the realizations this chunk must process
			realizations_in_chunk = realization_list[realizations_per_chunk*c:realizations_per_chunk*(c+1)]

			#Find in how many executables we should split the realizations in this chunk
			num_executables = len(realizations_in_chunk)

			#Find how many realizations per executable
			executables = list()
			cores = [job_settings.cores_per_simulation] * num_executables
			
			for e in range(num_executables):
			
				#Separate the cosmo_id,geometry_id,realization number
				cosmo_id,geometry_id,ic_number = realizations_in_chunk[e].split("|")

				#Get the corresponding SimulationXXX instances
				model = self.getModel(cosmo_id)

				nside,box_size = geometry_id.split("b")
				nside = int(nside)
				box_size = float(box_size) * model.Mpc_over_h
				collection = model.getCollection(box_size=box_size,nside=nside)

				r = collection.getRealization(int(ic_number.strip("ic")))

				#Check that the cores per simulation matches the number of files per snapshot
				with self.syshandler.open(os.path.join(r.home_subdir,"gadget2.p"),"rb") as settingsfile:
					gadget_settings = self.syshandler.pickleload(settingsfile)

				assert gadget_settings.NumFilesPerSnapshot==job_settings.cores_per_simulation,"In the current implementation of plane generation, the number of MPI tasks must be the same as the number of files per snapshot!"

				#Figure out the correct environment file
				if "environment_file" in kwargs.keys():
					environment_file = kwargs["environment_file"]
				else:
					environment_file = "environment.ini"

				#Figure out the correct configuration file
				if "plane_config_file" in kwargs.keys():
					config_file = kwargs["plane_config_file"]
				else:
					config_file = "lens.ini"

				executables.append(job_settings.path_to_executable + " " + """-e {0} -c {1} "{2}" """.format(environment_file,config_file,realization_list[realizations_per_chunk*c+e]))

			#Write the script
			script_filename = os.path.join(self.environment.home,"Jobs",job_settings.job_script_file)

			if not one_script:
				script_filename_split = script_filename.split(".")
				script_filename_split[-2] += "{0}".format(c+1)
				script_filename = ".".join(script_filename_split)

			#Override settings to make stdout and stderr go in the right places
			job_settings.redirect_stdout = os.path.join(self.environment.home,"Logs",job_settings.redirect_stdout)
			job_settings.redirect_stderr = os.path.join(self.environment.home,"Logs",job_settings.redirect_stderr)

			if (not one_script) or (not c):

				#Inform user where logs will be directed
				print("[+] Stdout will be directed to {0}".format(job_settings.redirect_stdout))
				print("[+] Stderr will be directed to {0}".format(job_settings.redirect_stderr))

				with self.syshandler.open(script_filename,"w") as scriptfile:
					scriptfile.write(job_handler.writePreamble(job_settings))
			
			with self.syshandler.open(script_filename,"a") as scriptfile:
				scriptfile.write(job_handler.writeExecution(executables,cores,job_settings))

			#Log to user and return
			if (not one_script) or (not c):
				print("[+] {0} written on {1}".format(script_filename,self.syshandler.name))

	############################################################################################################################################

	def writeRaySubmission(self,realization_list,job_settings,job_handler,chunks=1,**kwargs):
		
		"""
		Writes raytracing submission script

		:param realization_list: list of ics to generate in the form "cosmo_id|geometry_id"
		:type realization_list: list. of str.

		:param chunks: number of independent jobs in which to split the submission (one script per job will be written)
		:type chunks: int. 

		:param job_settings: settings for the job (resources, etc...)
		:type job_settings: JobSettings

		:param job_handler: handler of the cluster specific features (job scheduler, architecture, etc...)
		:type job_handler: JobHandler

		:param kwargs: keyword arguments accepted are "environment_file" to specify the environment settings for the current batch, "raytracing_config_file" to specify the lensing option for the ray tracing. Additionally you can set one_script=True to include all the executables sequentially in a single script
		:type kwargs: dict.

		"""

		#Type safety check
		assert isinstance(job_settings,JobSettings)
		assert isinstance(job_handler,JobHandler)
		assert len(realization_list)%chunks==0,"Perfect load balancing enforced, each job will process the same number of realizations!"

		#Check if we need to collapse everyting in one script
		if "one_script" in kwargs.keys():
			one_script = kwargs["one_script"]
		else:
			one_script = False

		#Create the dedicated Job and Logs directories if not existent already
		for d in [os.path.join(self.environment.home,"Jobs"),os.path.join(self.environment.home,"Logs")]:
			
			if not(self.syshandler.exists(d)):
				self.syshandler.mkdir(d)
				print("[+] {0} created on {1}".format(d,self.syshandler.name))

		#Split realizations between independent jobs
		realizations_per_chunk = len(realization_list)//chunks

		#Override job settings to make sure requested resources are enough
		job_settings.num_cores = job_settings.cores_per_simulation*realizations_per_chunk

		for c in range(chunks):

			#Find the realizations this chunk must process
			realizations_in_chunk = realization_list[realizations_per_chunk*c:realizations_per_chunk*(c+1)]

			#Find in how many executables we should split the realizations in this chunk
			num_executables = len(realizations_in_chunk)

			#Find how many realizations per executable
			executables = list()
			cores = [job_settings.cores_per_simulation] * num_executables
			
			for e in range(num_executables):
			
				#Separate the cosmo_id,geometry_id,realization number
				parts = realizations_in_chunk[e].split("|")

				#Figure out the correct environment file
				if "environment_file" in kwargs.keys():
					environment_file = kwargs["environment_file"]
				else:
					environment_file = "environment.ini"

				#Figure out the correct configuration file
				if "raytracing_config_file" in kwargs.keys():
					config_file = kwargs["raytracing_config_file"]
				else:
					config_file = "configuration.ini"

				#Make sure that there will be perfect load balancing at execution
				if len(parts)==2:

					try:
						raytracing_settings = MapSettings.read(config_file)
						if raytracing_settings.lens_map_realizations%job_settings.cores_per_simulation:
							raise ValueError("The number of map realizations must be a multiple of the number of cores per simulation!")

					except AssertionError:
						raytracing_settings = CatalogSettings.read(config_file)
						if raytracing_settings.lens_catalog_realizations%job_settings.cores_per_simulation:
							raise ValueError("The number of map realizations must be a multiple of the number of cores per simulation!")

				elif len(parts)==1:
					
					raytracing_settings = TelescopicMapSettings.read(config_file)
					if raytracing_settings.lens_map_realizations%job_settings.cores_per_simulation:
						raise ValueError("The number of map realizations must be a multiple of the number of cores per simulation!")

				else:
					raise ValueError("There are too many '|'' into your id: {0}".format(realizations_in_chunk[e]))

				executables.append(job_settings.path_to_executable + " " + """-e {0} -c {1} "{2}" """.format(environment_file,config_file,realization_list[realizations_per_chunk*c+e]))

			#Write the script
			script_filename = os.path.join(self.environment.home,"Jobs",job_settings.job_script_file)

			if not one_script:
				script_filename_split = script_filename.split(".")
				script_filename_split[-2] += "{0}".format(c+1)
				script_filename = ".".join(script_filename_split)

			#Override settings to make stdout and stderr go in the right places
			job_settings.redirect_stdout = os.path.join(self.environment.home,"Logs",job_settings.redirect_stdout)
			job_settings.redirect_stderr = os.path.join(self.environment.home,"Logs",job_settings.redirect_stderr)


			if (not one_script) or (not c):
			
				#Inform user where logs will be directed
				print("[+] Stdout will be directed to {0}".format(job_settings.redirect_stdout))
				print("[+] Stderr will be directed to {0}".format(job_settings.redirect_stderr))

				with self.syshandler.open(script_filename,"w") as scriptfile:
					scriptfile.write(job_handler.writePreamble(job_settings))

			with self.syshandler.open(script_filename,"a") as scriptfile:
				scriptfile.write(job_handler.writeExecution(executables,cores,job_settings))

			#Log to user and return
			if (not one_script) or (not c):
				print("[+] {0} written on {1}".format(script_filename,self.syshandler.name))

	############################################################################################################################################

	def writeSubmission(self,realization_list,job_settings,job_handler,job_executable="lenstools.custom",chunks=1,**kwargs):
		
		"""
		Writes a generic job submission script

		:param realization_list: list of ics to generate in the form "cosmo_id|geometry_id"
		:type realization_list: list. of str.

		:param chunks: number of independent jobs in which to split the submission (one script per job will be written)
		:type chunks: int. 

		:param job_settings: settings for the job (resources, etc...)
		:type job_settings: JobSettings

		:param job_handler: handler of the cluster specific features (job scheduler, architecture, etc...)
		:type job_handler: JobHandler

		:param job_executable: name of the executable that will be run
		:type job_executable: str.

		:param kwargs: keyword arguments accepted are "environment_file" to specify the environment settings for the current batch, "config_file" to specify the job specific configuration file. Additionally you can set one_script=True to include all the executables sequentially in a single script
		:type kwargs: dict.

		"""

		#Type safety check
		assert isinstance(job_settings,JobSettings)
		assert isinstance(job_handler,JobHandler)
		assert len(realization_list)%chunks==0,"Perfect load balancing enforced, each job will process the same number of realizations!"

		#Check if we need to collapse everyting in one script
		if "one_script" in kwargs.keys():
			one_script = kwargs["one_script"]
		else:
			one_script = False

		#Create the dedicated Job and Logs directories if not existent already
		for d in [os.path.join(self.environment.home,"Jobs"),os.path.join(self.environment.home,"Logs")]:
			
			if not(self.syshandler.exists(d)):
				self.syshandler.mkdir(d)
				print("[+] {0} created on {1}".format(d,self.syshandler.name))

		#Split realizations between independent jobs
		realizations_per_chunk = len(realization_list)//chunks

		#Override job settings to make sure requested resources are enough
		job_settings.num_cores = job_settings.cores_per_simulation*realizations_per_chunk

		for c in range(chunks):

			#Find the realizations this chunk must process
			realizations_in_chunk = realization_list[realizations_per_chunk*c:realizations_per_chunk*(c+1)]

			#Find in how many executables we should split the realizations in this chunk
			num_executables = len(realizations_in_chunk)

			#Find how many realizations per executable
			executables = list()
			cores = [job_settings.cores_per_simulation] * num_executables
			
			for e in range(num_executables):

				#Figure out the correct environment file
				if "environment_file" in kwargs.keys():
					environment_file = kwargs["environment_file"]
				else:
					environment_file = "environment.ini"

				#Figure out the correct configuration file
				if "config_file" in kwargs.keys():
					config_file = kwargs["config_file"]
				else:
					config_file = "configuration.ini"

				executables.append(job_executable + " " + """-e {0} -c {1} "{2}" """.format(environment_file,config_file,realization_list[realizations_per_chunk*c+e]))

				#Allow for extra arguments in the executable
				if "extra_args" in kwargs.keys() and (kwargs["extra_args"] is not None):
					executables[-1] += kwargs["extra_args"]

			#Write the script
			script_filename = os.path.join(self.environment.home,"Jobs",job_settings.job_script_file)

			if not one_script:
				script_filename_split = script_filename.split(".")
				script_filename_split[-2] += "{0}".format(c+1)
				script_filename = ".".join(script_filename_split)

			#Override settings to make stdout and stderr go in the right places
			job_settings.redirect_stdout = os.path.join(self.environment.home,"Logs",job_settings.redirect_stdout)
			job_settings.redirect_stderr = os.path.join(self.environment.home,"Logs",job_settings.redirect_stderr)


			if (not one_script) or (not c):
			
				#Inform user where logs will be directed
				print("[+] Stdout will be directed to {0}".format(job_settings.redirect_stdout))
				print("[+] Stderr will be directed to {0}".format(job_settings.redirect_stderr))

				with self.syshandler.open(script_filename,"w") as scriptfile:
					scriptfile.write(job_handler.writePreamble(job_settings))

			with self.syshandler.open(script_filename,"a") as scriptfile:
				scriptfile.write(job_handler.writeExecution(executables,cores,job_settings))

			#Log to user and return
			if (not one_script) or (not c):
				print("[+] {0} written on {1}".format(script_filename,self.syshandler.name))	


##############################################
##############TreeNode class##################
##############################################

class TreeNode(object):

	__metaclass__ = ABCMeta

	@abstractmethod
	def __init__(self,*args):
		pass

	#######
	#Paths#
	#######

	@property
	def home(self):
		return self.home_subdir

	@property
	def storage(self):
		return self.syshandler.map(self.storage_subdir)

	@property
	def infofile(self):
		return os.path.join(self.environment.home,self.environment.json_tree_file)

	def path(self,filename,where="storage_subdir"):

		"""
		Returns the complete path to the resource corresponding to filename; returns None if the resource does not exist

		:param filename: name of the resource
		:type filename: str.

		:param where: where to look for the resource, can be "home_subdir" or "storage_subdir"
		:type where: str.

		:returns: full path to the resource
		:rtype: str.

		"""

		return self.syshandler.map(os.path.join(getattr(self,where),filename))

	def ls(self,glob="*",where="storage_subdir"):

		"""
		Returns the list of files present either in the storage or home portion of the simulation batch

		:param glob: glob string to filter file types
		:type glob: str. 

		:param where: specifies if to look into the storage or home part of the simulation batch
		:type where: str.

		:returns: list of files
		:rtype: list.

		"""

		search_path = getattr(self,where)
		return [os.path.relpath(f,search_path) for f in self.syshandler.glob(os.path.join(search_path,glob))]

	################################################################################################################################

	def mkdir(self,directory):

		"""
		Create a sub-directory inside the current instance home and storage paths

		:param directory: name of the directory to create
		:type directory: str.

		"""

		for d in [self.home_subdir,self.storage_subdir]:

			dir_to_make = os.path.join(d,directory)
			
			if not self.syshandler.exists(dir_to_make):
				self.syshandler.mkdir(dir_to_make)
				print("[+] {0} created on {1}".format(dir_to_make,self.syshandler.name))

	################################################################################################################################

	def mkfifo(self,names,method=None,where="storage_subdir"):

		"""
		Makes a list of named pipes

		:param names: names of the named pipes
		:type names: list.

		:param where: where to put the named pipes, default is "storage_subdir"
		:type where: str.

		"""

		if method is None:
			method = lambda n:os.mkfifo(n)

		for name in names:
			pipe_name = self.syshandler.map(os.path.join(getattr(self,where),name))
			method(pipe_name)
			print("[+] Created named pipe {0}".format(pipe_name)) 

	################################################################################################################################

	def execute(self,filename,callback=None,where="storage_subdir",**kwargs):

		"""
		Calls a user defined function on the resource pointed to by filename; if None is provided, returns the full path to the map

		:param filename: name of the file on which to call the callback
		:type filename: str.

		:param callback: user defined function that takes filename as first argument
		:type callback: callable

		:param where: where to look for the resource, can be "home_subdir" or "storage_subdir"
		:type where: str.

		:param kwargs: key word arguments to be passed to callback
		:type kwargs: dict.

		:returns: the result of callback

		"""

		full_path = self.path(filename,where)
		if full_path is None:
			raise IOError("{0} does not exist!".format(filename))

		if callback is None:
			return full_path

		#Call the function
		return callback(full_path,**kwargs)

#####################################################
##############SimulationModel class##################
#####################################################

class SimulationModel(TreeNode):

	"""
	Class handler of a weak lensing simulation model, defined by a set of cosmological parameters

	"""

	@property
	def cosmo_id(self):
		base = "{0}{1:."+str(self.environment.cosmo_id_digits)+"f}"
		return "_".join([ base.format(p,getattr(self.cosmology,self.environment.name2attr[p])) for p in self.parameters if (hasattr(self.cosmology,self.environment.name2attr[p]) and getattr(self.cosmology,self.environment.name2attr[p]) is not None)])

	@property
	def info(self):
		with InfoDict(self.environment,syshandler=self.syshandler) as infodict:
			return infodict.dictionary

	def __init__(self,cosmology=configuration.CosmoDefault,environment=None,parameters=["Om","Ol","w","ns","si"],**kwargs):

		"""
		Set the base for the simulation

		:param cosmology: cosmological model to simulate
		:type cosmology: LensToolsCosmology

		:param environment: environment settings of the current machine
		:type environment: EnvironmentSettings

		:param parameters: cosmological parameters to keep track of
		:type parameters: list.

		"""

		#Safety checks
		assert isinstance(cosmology,LensToolsCosmology)

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

		#Create directories accordingly
		self.home_subdir = os.path.join(self.environment.home,self.cosmo_id)
		self.storage_subdir = os.path.join(self.environment.storage,self.cosmo_id)

		for key in kwargs.keys():
			setattr(self,key,kwargs[key])
	

	def __repr__(self):

		representation_parameters = []
		for p in self.parameters:
			representation_parameters.append(("{0}={1:."+str(self.environment.cosmo_id_digits)+"f}").format(p,getattr(self.cosmology,self.environment.name2attr[p])))

		return "<"+ " , ".join(representation_parameters) + ">"


	#Convenient resource retrieval
	def __getitem__(self,s):
		return _retrieve(self,s)

	################################################################################################################################

	def plan(self,map_size,max_redshift,nlenses,mass_resolution=1.0e10*u.Msun):

		"""
		This method is useful when planning a simulation batch: given an expected size of the simulated maps, a maximum redshift for the ray tracing and the number of lenses along the line of sight, it computes the size of the box that must be employed, the spacing between the lenses and the redshifts at which the lenses need to be placed

		:param map_size: angular size of the simulated maps
		:type map_size: quantity

		:param max_redshift: maximum redshift of the maps
		:type max_redshift: float.

		:param nlenses: number of lenses along the line of sight
		:type nlenses: int.

		:param mass_resolution: mass of one particle in the N--body simulation
		:type mass_resoultion: quantity

		:returns: suggested simulations specifications
		:rtype: dict.

		"""

		#Comoving distance to max_redshift
		distance = self.cosmology.comoving_distance(max_redshift)

		#First compute the box size
		box_size = map_size.to(u.rad).value * distance

		#Now compute the total mass in the box
		total_mass = self.cosmology.critical_density0 * self.cosmology.Om0 * (box_size**3)
		nside = ((total_mass / mass_resolution).decompose().value)**(1/3)

		#Next compute the lens plane thickness
		thickness = distance / nlenses

		#Compute the redshifts at which the lenses need to be put in
		lens_distance = np.linspace(thickness.value,distance.value,nlenses) * distance.unit
		lens_redshift = np.zeros(len(lens_distance))

		for n,d in enumerate(lens_distance):
			lens_redshift[n] = z_at_value(self.cosmology.comoving_distance,d)

		#Build the dictionary
		specs = dict()
		specs["box_size"] = box_size.to(self.Mpc_over_h)
		specs["nside"] = int(nside)
		specs["plane_thickness"] = thickness.to(self.Mpc_over_h)
		specs["output_scale_factor"] = np.sort(1/(1+lens_redshift))

		#Return
		return specs

	################################################################################################################################

	def newCollection(self,box_size=240.0*u.Mpc,nside=512):

		"""
		Instantiate new simulation with the specified settings

		"""

		newSimulation = SimulationCollection(self.cosmology,self.environment,self.parameters,box_size,nside,syshandler=self.syshandler)

		#Make the corresponding directory if not already present
		for d in [newSimulation.home_subdir,newSimulation.storage_subdir]:
			if not self.syshandler.exists(d):

				self.syshandler.mkdir(d)
				print("[+] {0} created on {1}".format(d,self.syshandler.name))

				#Update dictionary if the batch is indicized
				if self.syshandler.exists(self.infofile):
					self.info[newSimulation.cosmo_id][newSimulation.geometry_id] = dict()
					self.info[newSimulation.cosmo_id][newSimulation.geometry_id]["nbody"] = dict()
					self.info[newSimulation.cosmo_id][newSimulation.geometry_id]["map_sets"] = dict()
					self.info[newSimulation.cosmo_id][newSimulation.geometry_id]["catalogs"] = dict()

			else:
				print("[-] Collection {0} already exists!".format(os.path.join(newSimulation.cosmo_id,newSimulation.geometry_id)))
				return newSimulation
				

		#Keep track of the fact we created a new collection
		with self.syshandler.open(os.path.join(self.environment.home,"collections.txt"),"a") as logfile:
			logfile.write("{0}|{1}\n".format(self.cosmo_id,newSimulation.geometry_id))

		return newSimulation

	################################################################################################################################

	def getCollection(self,box_size,nside=None):

		#Allow to pass a geometry_id as first argument
		if isinstance(box_size,str):
			try:
				parts = box_size.split("b")
				nside = int(parts[0])
				box_size = float(parts[1]) * self.Mpc_over_h
			except (IndexError,ValueError):
				return None

		assert nside is not None,"if you did not specify the second argument, it means the first should be in the form 'xxxbyyy'"

		#See if the collection exists
		collection = SimulationCollection(self.cosmology,self.environment,self.parameters,box_size,nside,syshandler=self.syshandler)
		if not(self.syshandler.exists(collection.storage_subdir) and self.syshandler.exists(collection.home_subdir)):
			return None

		#If it exists, return the SimulationCollection instance
		return collection

	################################################################################################################################


	@property
	def collections(self):

		"""
		Lists all the available collections for a model

		:returns: list.
		:rtype: SimulationCollection

		"""

		if self.syshandler.exists(self.infofile):
			collection_names = self.info[self.cosmo_id].keys()
		else:
			collection_names = [ os.path.basename(d) for d in self.syshandler.glob(os.path.join(self.home_subdir,"*")) ]

		collection_list = list()

		for name in collection_names:
			collection = self.getCollection(name)
			if collection is not None:
				collection_list.append(collection)

		return collection_list


	################################################
	############Telescopic map sets#################
	################################################

	def newTelescopicMapSet(self,collections,redshifts,settings):

		"""
		Create a new telescopic map set with the provided settings

		:param collections: list of the SimulationCollection instances that participate in the telescopic map set
		:type collections: list.

		:param redshifts: redshifts that mark the transition from one collection to the other during raytracing
		:type redshifts: array.

		:param settings: settings for the new map set
		:type settings: TelescopicMapSettings

		:rtype: SimulationMaps

		"""

		#Safety check
		assert isinstance(settings,TelescopicMapSettings)
		assert type(settings.plane_set)==tuple
		assert type(settings.mix_nbody_realizations)==tuple
		assert type(settings.mix_cut_points)==tuple
		assert type(settings.mix_normals)==tuple

		#Adjust according to number of collections (one specification for each different collection to use)
		for attr_name in ["plane_set","mix_nbody_realizations","mix_cut_points","mix_normals"]:
			
			attr = getattr(settings,attr_name)

			if len(attr)!=len(collections):
				if len(attr)==1:
					setattr(settings,attr_name,attr*len(collections))
				else:
					raise ValueError("You need to specify a {0} for each collection!".format(attr_name)) 


		#Instantiate SimulationMaps object
		map_set = SimulationTelescopicMaps(self.cosmology,self.environment,self.parameters,collections,redshifts,settings,syshandler=self.syshandler)

		#Create dedicated directories
		for d in [ map_set.home_subdir,map_set.storage_subdir ]:
			if not self.syshandler.exists(d):
				self.syshandler.mkdir(d)
				print("[+] {0} created on {1}".format(d,self.syshandler.name))

		#Save a picked copy of the settings to use for future reference
		with self.syshandler.open(os.path.join(map_set.home_subdir,"settings.p"),"wb") as settingsfile:
			self.syshandler.pickledump(settings,settingsfile)

		#Append the name of the map batch to a summary file
		with self.syshandler.open(os.path.join(self.home_subdir,"telescopic_map_sets.txt"),"a") as setsfile:
			setsfile.write("{0},{1},{2}\n".format(map_set.name,"-".join([c.geometry_id for c in map_set.mapcollections]),"-".join(["{0:.3f}".format(z) for z in map_set.redshifts])))

		#Return to user
		return map_set

	################################################################################################################################

	def getTelescopicMapSet(self,setname):

		"""
		Return an instance of the telescopic map set named "name"

		:param setname: name of the map set
		:type setname: str.

		:rtype: SimulationTelescopicMaps

		"""

		with self.syshandler.open(os.path.join(self.home_subdir,"telescopic_map_sets.txt"),"r") as setsfile:
			
			while True:
				
				line = setsfile.readline().strip("\n")
				if line=="":
					return None

				name,colnames,rednames = line.rstrip("\n").split(",")

				if name!=setname:
					continue

				#Pickle the settings
				with self.syshandler.open(os.path.join(self.home_subdir,name,"settings.p"),"rb") as settingsfile:
					settings = self.syshandler.pickleload(settingsfile)

				#Parse the collections
				collections = [ self.getCollection(colname) for colname in colnames.split("-") ]

				#Parse the redshifts
				redshifts = np.array([ float(z) for z in rednames.split("-") ])

				#Return
				return SimulationTelescopicMaps(self.cosmology,self.environment,self.parameters,collections,redshifts,settings,syshandler=self.syshandler)


	################################################################################################################################

	@property
	def telescopicmapsets(self):

		"""
		Lists all the available telescopic map sets for a model

		:returns: list.
		:rtype: SimulationTelescopicMaps

		"""

		map_sets = list()

		try:
			with self.syshandler.open(os.path.join(self.home_subdir,"telescopic_map_sets.txt"),"r") as setsfile:
				lines = setsfile.readlines()
		except IOError:
			return map_sets

		for line in lines:

			name,colnames,rednames = line.rstrip("\n").split(",")

			#Pickle the settings
			with self.syshandler.open(os.path.join(self.home_subdir,name,"settings.p"),"rb") as settingsfile:
				settings = self.syshandler.pickleload(settingsfile)

			#Parse the collections
			collections = [ self.getCollection(colname) for colname in colnames.split("-") ]

			#Parse the redshifts
			redshifts = np.array([ float(z) for z in rednames.split("-") ])

			#Append the telescopic map instance to map_sets
			map_sets.append(SimulationTelescopicMaps(self.cosmology,self.environment,self.parameters,collections,redshifts,settings,syshandler=self.syshandler))

		#Return to user
		return map_sets


##########################################################
##############SimulationCollection class##################
##########################################################

class SimulationCollection(SimulationModel):


	"""
	Class handler of a collection of simulations that share model parameters

	"""

	def __init__(self,cosmology,environment,parameters,box_size,nside,**kwargs):

		super(SimulationCollection,self).__init__(cosmology,environment,parameters,**kwargs)
		self.box_size = box_size
		self.nside = nside

		#Build the geometry_id
		self.geometry_id = "{0}b{1}".format(nside,int(box_size.to(self.Mpc_over_h).value))

		#Build the directory names
		self.home_subdir = os.path.join(self.environment.home,self.cosmo_id,self.geometry_id)
		self.storage_subdir = os.path.join(self.environment.storage,self.cosmo_id,self.geometry_id)


	def __repr__(self):

		return super(SimulationCollection,self).__repr__() + " | box={0},nside={1}".format(self.box_size,self.nside)


	################################################################################################################################

	#Can't call these methods anymore
	def newCollection(self):
		raise NotImplementedError

	def getCollection(self):
		raise NotImplementedError

	@property
	def collections(self):
		raise NotImplementedError

	def newTelescopicMapSet(self):
		raise NotImplementedError

	def getTelescopicMapSet(self):
		raise NotImplementedError

	@property
	def telescopicmapsets(self):
		raise NotImplementedError		

	################################################################################################################################

	@property
	def resolution(self):

		"""
		Computes the mass resolution (mass of one particle) of the simulation collection

		"""

		dm_density = self.cosmology.critical_density0 * self.cosmology.Om0
		return (dm_density*(self.box_size**3) / (self.nside**3)).to(u.Msun)

	################################################################################################################################

	def tile(self,n):

		"""
		Computes the location of the slab centers and the thickness of the slabs necessary to cut the simulation box in n equal pieces

		:param n: the number of pieces to cut the box into
		:type n: int.

		:return: thickness,cut_points
		:rtype: tuple. 

		"""

		#Compute thickness
		thickness = self.box_size / n

		#Compute cut points
		cut_points = np.zeros(n) * self.Mpc_over_h

		for p in range(n):
			cut_points[p] = thickness*(0.5 + p)

		#Return in physical units
		return thickness.to(u.Mpc),cut_points.to(u.Mpc)

	################################################################################################################################

	def newRealization(self,seed=0,**kwargs):
		
		#Check if there are already generated initial conditions in there
		ics_present = self.syshandler.glob(os.path.join(self.storage_subdir,"ic*"))
		new_ic_index = len(ics_present) + 1

		#Generate the new initial condition
		newIC = SimulationIC(self.cosmology,self.environment,self.parameters,self.box_size,self.nside,new_ic_index,seed,syshandler=self.syshandler,**kwargs)

		#Make dedicated directories for new initial condition,ics and snapshots
		for d in [newIC.home_subdir,newIC.storage_subdir,newIC.ics_subdir,newIC.snapshot_subdir]:

			if not self.syshandler.exists(d):
				self.syshandler.mkdir(d)
				print("[+] {0} created on {1}".format(d,self.syshandler.name))

				#Update dictionary if the batch is indicized
				if self.syshandler.exists(self.infofile):
					self.info[newIC.cosmo_id][newIC.geometry_id]["nbody"][str(newIC.ic_index)] = dict()
					self.info[newIC.cosmo_id][newIC.geometry_id]["nbody"][str(newIC.ic_index)]["plane_sets"] = dict()

		#Make new file with the number of the seed, write the ICFileBase, SnapshotFileBase to the seed file
		with self.syshandler.open(os.path.join(newIC.home_subdir,"seed"+str(seed)),"w") as seedfile:
			contents = { "ICFileBase":newIC.ICFileBase, "SnapshotFileBase":newIC.SnapshotFileBase }
			seedfile.write(json.dumps(contents))

		#Keep track of the fact that we created a new nbody realization
		with self.syshandler.open(os.path.join(self.environment.home,"realizations.txt"),"a") as logfile:
			logfile.write("{0}|{1}|ic{2}\n".format(self.cosmo_id,self.geometry_id,new_ic_index))

		return newIC

	################################################################################################################################

	def getRealization(self,n):

		#Check if this particular realization exists
		newIC = SimulationIC(self.cosmology,self.environment,self.parameters,self.box_size,self.nside,n,seed=0,syshandler=self.syshandler)
		if not(self.syshandler.exists(newIC.home_subdir)) or not(self.syshandler.exists(newIC.storage_subdir)):
			return None

		#Read the seed number
		seed_filename = self.syshandler.glob(os.path.join(newIC.home_subdir,"seed*"))[0]
		newIC.seed = int(os.path.basename(seed_filename).strip("seed"))

		#Read ICFileBase,SnapshotFileBase from seed file
		with self.syshandler.open(seed_filename,"r") as fp:

			try:
				contents = json.loads(fp.read())
				newIC.ICFileBase = contents["ICFileBase"]
				newIC.SnapshotFileBase = contents["SnapshotFileBase"]
			except ValueError:
				pass

		#Return the SimulationIC instance
		return newIC

	################################################################################################################################

	@property
	def realizations(self):

		"""
		List the available realizations (or independent simulations) for the current collection

		:returns: list.
		:rtype: SimulationIC

		"""

		if self.syshandler.exists(self.infofile):
			ic_numbers = [ int(n) for n in self.info[self.cosmo_id][self.geometry_id]["nbody"] ]
		else:
			ic_numbers = [ int(os.path.basename(d).strip("ic")) for d in self.syshandler.glob(os.path.join(self.home_subdir,"ic*")) ]
		
		#Get realizations
		ic_list = [ self.getRealization(ic) for ic in sorted(ic_numbers) ]

		#Return to user
		return ic_list


	################################################################################################################################

	def newMapSet(self,settings):

		"""
		Create a new map set with the provided settings

		:param settings: settings for the new map set
		:type settings: MapSettings

		:rtype: SimulationMaps

		"""

		#Safety check
		assert isinstance(settings,MapSettings)

		#Instantiate SimulationMaps object
		map_set = SimulationMaps(self.cosmology,self.environment,self.parameters,self.box_size,self.nside,settings,syshandler=self.syshandler)

		#Create dedicated directories
		for d in [ map_set.home_subdir,map_set.storage_subdir ]:
			if not self.syshandler.exists(d):
				self.syshandler.mkdir(d)
				print("[+] {0} created on {1}".format(d,self.syshandler.name))

				#Update dictionary if the batch is indicized
				if self.syshandler.exists(self.infofile):
					self.info[map_set.cosmo_id][map_set.geometry_id]["map_sets"][settings.directory_name] = dict()

		#Save a picked copy of the settings to use for future reference
		with self.syshandler.open(os.path.join(map_set.home_subdir,"settings.p"),"wb") as settingsfile:
			self.syshandler.pickledump(settings,settingsfile)

		#Append the name of the map batch to a summary file
		if not self.syshandler.exists(self.infofile):
			with self.syshandler.open(os.path.join(self.home_subdir,"sets.txt"),"a") as setsfile:
				setsfile.write("{0}\n".format(settings.directory_name))

		#Return to user
		return map_set

	#################################################################################################################################

	def getMapSet(self,setname):

		"""
		Get access to a pre-existing map set

		:param setname: name of the map set to access
		:type setname: str.

		"""

		#Check if the map set exists
		if (not self.syshandler.exists(os.path.join(self.storage_subdir,setname))) or (not self.syshandler.exists(os.path.join(self.home_subdir,setname))):
			return None

		#Read the settings from the pickled file
		with self.syshandler.open(os.path.join(self.home_subdir,setname,"settings.p"),"rb") as settingsfile:
			settings = self.syshandler.pickleload(settingsfile) 

		#Return to user
		return SimulationMaps(self.cosmology,self.environment,self.parameters,self.box_size,self.nside,settings,syshandler=self.syshandler)

	#################################################################################################################################

	@property
	def mapsets(self):

		"""
		Build a list with all the available map sets

		:returns: list of SimulationMaps

		"""

		map_sets = list()
		
		try:

			if self.syshandler.exists(self.infofile):
				for set_name in self.info[self.cosmo_id][self.geometry_id]["map_sets"].keys():
					map_sets.append(self.getMapSet(set_name))
			else:
				with self.syshandler.open(os.path.join(self.home_subdir,"sets.txt"),"r") as setsfile:
					for set_name in setsfile.readlines():
						map_sets.append(self.getMapSet(set_name.strip("\n")))

		except IOError:
			pass
		finally:
			return map_sets

	#################################################################################################################################

	def newCatalog(self,settings):

		"""
		Create a new simulated catalog with the provided settings

		:param settings: settings for the new simulated catalog
		:type settings: CatalogSettings

		:rtype: SimulationCatalog

		"""

		#Safety check
		assert isinstance(settings,CatalogSettings)

		#Instantiate SimulationMaps object
		catalog = SimulationCatalog(self.cosmology,self.environment,self.parameters,self.box_size,self.nside,settings,syshandler=self.syshandler)

		#Create dedicated directories
		for d in [ catalog.home_subdir,catalog.storage_subdir ]:
			if not self.syshandler.exists(d):
				self.syshandler.mkdir(d)
				print("[+] {0} created on {1}".format(d,self.syshandler.name))

				#Update dictionary if the batch is indicized
				if self.syshandler.exists(self.infofile):
					self.info[catalog.cosmo_id][catalog.geometry_id]["catalogs"][settings.directory_name] = dict()

		#Save a picked copy of the settings to use for future reference
		with self.syshandler.open(os.path.join(catalog.home_subdir,"settings.p"),"wb") as settingsfile:
			self.syshandler.pickledump(settings,settingsfile)

		#Append the name of the map batch to a summary file
		if not self.syshandler.exists(self.infofile):
			with self.syshandler.open(os.path.join(self.home_subdir,"catalogs.txt"),"a") as setsfile:
				setsfile.write("{0}\n".format(settings.directory_name))

		#Return to user
		return catalog

	#################################################################################################################################

	def getCatalog(self,catalog_name):

		"""
		Get access to a pre-existing catalog

		:param catalog_name: name of the catalog to access
		:type catalog_name: str.

		"""

		#Check if the map set exists
		if (not self.syshandler.exists(os.path.join(self.storage_subdir,catalog_name))) or (not self.syshandler.exists(os.path.join(self.home_subdir,catalog_name))):
			return None

		#Read the settings from the pickled file
		with self.syshandler.open(os.path.join(self.home_subdir,catalog_name,"settings.p"),"rb") as settingsfile:
			settings = self.syshandler.pickleload(settingsfile) 

		#Return to user
		return SimulationCatalog(self.cosmology,self.environment,self.parameters,self.box_size,self.nside,settings,syshandler=self.syshandler)


	#################################################################################################################################

	@property
	def catalogs(self):

		"""
		Builds a list with all the available catalogs

		:returns: list of SimulationCatalog

		"""

		catalogs = list()
		
		try:
			
			if self.syshandler.exists(self.infofile):
				for catalog_name in self.info[self.cosmo_id][self.geometry_id]["catalogs"].keys():
					catalogs.append(self.getCatalog(catalog_name))
			else:
				with self.syshandler.open(os.path.join(self.home_subdir,"catalogs.txt"),"r") as setsfile:
					for catalog_name in setsfile.readlines():
						catalogs.append(self.getCatalog(catalog_name.strip("\n")))
		
		except IOError:
			pass
		finally:
			return catalogs


	###################################################################
	##################Useful methods for CAMB I/O######################
	###################################################################

	def writeCAMB(self,z,settings,fname="camb.param",output_root="camb"):

		"""
		Generates the parameter file that CAMB needs to read to evolve the current initial condition in time

		:param settings: CAMB tunable settings
		:type settings: CAMBSettings

		:param z: redshifts at which CAMB needs to compute the matter power spectrum
		:type z: array.

		:param fname: name of the parameter file to write
		:type fname: str.

		:param output_root: output root of camb products
		:type output_root: str.

		"""

		#Safety type check
		assert isinstance(settings,CAMBSettings)
		if type(z)==np.float:
			z = np.array([z])

		#Write the parameter file
		camb_filename = os.path.join(self.home_subdir,fname) 
		with self.syshandler.open(camb_filename,"w") as paramfile:
			paramfile.write(settings.write(output_root=os.path.join(self.home_subdir,output_root),cosmology=self.cosmology,redshifts=z))

		print("[+] {0} written on {1}".format(camb_filename,self.syshandler.name))

		#Save a pickled copy of the settings
		with self.syshandler.open(os.path.join(self.home_subdir,"camb.p"),"wb") as settingsfile:
			self.syshandler.pickledump(settings,settingsfile)

	def loadTransferFunction(self,name_root="camb_matterpower"):

		"""
		Loads in the CDM transfer function calculated with CAMB (k,delta(k)/k^2) for a unit super-horizon perturbation

		:param name_root: root of the file names that contain the transfer function
		:type name_root: str.

		:rtype: :py:class:`~lenstools.simulations.camb.CAMBTransferFromPower`

		"""

		#Look in the home subdirectory for files camb_transferfunc_zxxxxxxxx.dat
		tfr_files = self.ls(name_root+"_z*.dat","home")
		if not len(tfr_files):
			raise ValueError("Transfer functions have not been computed yet!")

		#Look for the wavenumbers in the first file
		k = np.loadtxt(self.path(tfr_files[0],"home"),usecols=(0,))*self.cosmology.h*(self.Mpc_over_h**-1)
		tfr = CAMBTransferFromPower(k)

		#Add transfer function information from each redshift
		for f in tfr_files:
			z, = re.search(r"z([0-9.]+)",f).groups()
			transfer_values = np.loadtxt(self.path(f,"home"),usecols=(1,))
			tfr.add(float(z.rstrip(".")),transfer_values)

		#Return to user
		return tfr

	def camb2ngenic(self,z,input_root="camb"):

		"""
		Read CAMB power spectrum file and convert it in a N-GenIC readable format

		:param z: redshift of the matter power spectrum file to convert
		:type z: float.

		:param input_root: name root of camb products
		:type input_root: str.

		"""

		camb_ps_file = os.path.join(self.environment.home,self.cosmo_id,self.geometry_id,"{0}_matterpower_z{1:.6f}.dat".format(input_root,z))
		if not(self.syshandler.exists(camb_ps_file)):
			raise IOError("CAMB matter power spectrum file {0} does not exist yet!!".format(camb_ps_file))

		k,Pk = np.loadtxt(camb_ps_file,unpack=True)
		lgk,lgP = _camb2ngenic(k,Pk)

		ngenic_ps_file = os.path.join(self.environment.home,self.cosmo_id,self.geometry_id,"ngenic_matterpower_z{0:.6f}.txt".format(z))
		np.savetxt(ngenic_ps_file,np.array([lgk,lgP]).T)

		print("[+] CAMB matter power spectrum at {0} converted into N-GenIC readable format at {1}".format(camb_ps_file,ngenic_ps_file))

################################################################################################################################


##########################################################
##############SimulationIC class##########################
##########################################################

class SimulationIC(SimulationCollection):

	"""
	Class handler of a simulation with a defined initial condition

	"""

	def __init__(self,cosmology,environment,parameters,box_size,nside,ic_index,seed,ICFileBase="ics",SnapshotFileBase="snapshot",**kwargs):

		super(SimulationIC,self).__init__(cosmology,environment,parameters,box_size,nside,**kwargs)

		#Save random seed information as attribute
		self.ic_index = ic_index
		self.seed = seed

		#Save storage sub-directory name
		self.home_subdir = os.path.join(self.home_subdir,"ic{0}".format(ic_index))
		self.storage_subdir = os.path.join(self.storage_subdir,"ic{0}".format(ic_index))

		#Save also snapshots and initial conditions dedicated sub-directories names
		self.ics_subdir = os.path.join(self.storage_subdir,"ics")
		self.snapshot_subdir = os.path.join(self.storage_subdir,"snapshots")

		#Useful to keep track of
		self.ICFileBase = ICFileBase
		self.SnapshotFileBase = SnapshotFileBase

		#Try to load in the simulation settings, if any are present
		try:
			with self.syshandler.open(os.path.join(self.home_subdir,"ngenic.p"),"rb") as settingsfile:
				self.ngenic_settings = self.syshandler.pickleload(settingsfile)
		except IOError:
			pass

		try:
			with self.syshandler.open(os.path.join(self.home_subdir,"gadget2.p"),"rb") as settingsfile:
				self.gadget_settings = self.syshandler.pickleload(settingsfile)
		except IOError:
			pass

	def __repr__(self):

		#Check if snapshots and/or initial conditions are present
		ics_on_disk = self.syshandler.glob(os.path.join(self.ics_subdir,self.ICFileBase+"*"))
		snap_on_disk = self.syshandler.glob(os.path.join(self.snapshot_subdir,self.SnapshotFileBase+"*"))

		return super(SimulationIC,self).__repr__() + " | ic={0},seed={1} | IC files on disk: {2} | Snapshot files on disk: {3}".format(self.ic_index,self.seed,len(ics_on_disk),len(snap_on_disk))

	#######
	#Paths#
	#######

	@property
	def ics(self):
		return self.syshandler.map(self.ics_subdir)

	@property
	def snapshots(self):
		return self.syshandler.map(self.snapshot_subdir)

	###########################################################################

	#Can't call these methods anymore
	def newRealization(self,seed):
		raise NotImplementedError

	@property
	def realizations(self):
		raise NotImplementedError

	def getRealization(self,ic):
		raise NotImplementedError

	def newMapSet(self,mp):
		raise NotImplementedError

	@property
	def mapsets(self):
		raise NotImplementedError

	def getMapSet(self,mp):
		raise NotImplementedError

	def newCatalog(self,ct):
		raise NotImplementedError

	@property
	def catalogs(self):
		raise NotImplementedError

	def getCatalog(self,ct):
		raise NotImplementedError

	def camb2ngenic(self,*args):
		raise NotImplementedError

	def writeCAMB(self,*args):
		raise NotImplementedError	

	####################################################################################################################################

	def snapshotPath(self,n,sub=0):

		"""
		Returns the complete path to a Gadget2 snapshot saved on disk; returns None if the snapshot does not exist

		:param n: number of the snapshot
		:type n: int.

		:param sub: number of the sub--file that contributes to the snapshot (if the snapshot has been split in multiple files) 
		:type sub: int.

		:returns: full path to snapshot
		:rtype: str.

		"""

		full_path = os.path.join(self.storage_subdir,"snapshots",self.gadget_settings.SnapshotFileBase+"_{0:03d}".format(n))

		if sub is None:
			return full_path

		if self.gadget_settings.NumFilesPerSnapshot > 1:
			full_path += ".{0}".format(sub)

		if not(self.syshandler.exists(full_path)):
			return None

		return full_path

	####################################################################################################################################

	def pipe_snapshots(self,nsnap=None,method=None):

		"""
		Create named pipes for nbody snapshots (to pipe them directly into the lens planes generation)

		:param nsnap: number of the snapshots to pipe
		:type nsnap: list.

		:param method: method to call on the filenames to make the named pipes
		:type method: callable

		"""

		files_per_snapshot = self.gadget_settings.NumFilesPerSnapshot
		if nsnap is None:
			nsnap = range(len(self.gadget_settings.OutputScaleFactor))

		if files_per_snapshot>1:
			snap_names = [ self.SnapshotFileBase+"_{0:03d}.{1}".format(i,j) for i,j in itertools.product(nsnap,range(files_per_snapshot)) ]
		else:
			snap_names = [ self.SnapshotFileBase+"_{0:03d}".format(n) for n in nsnap ] 

		self.mkfifo(snap_names,method=method,where="snapshot_subdir")


	####################################################################################################################################

	def newPlaneSet(self,settings):

		#Safety check
		assert isinstance(settings,PlaneSettings) or isinstance(settings,PlaneLightConeSettings)

		#Instantiate a SimulationPlanes object
		new_plane_set = SimulationPlanes(self.cosmology,self.environment,self.parameters,self.box_size,self.nside,self.ic_index,self.seed,self.ICFileBase,self.SnapshotFileBase,settings,syshandler=self.syshandler)

		#Create the dedicated directory if not present already
		for d in [new_plane_set.home_subdir,new_plane_set.storage_subdir]:
			if not(self.syshandler.exists(d)):
				self.syshandler.mkdir(d)
				print("[+] {0} created on {1}".format(d,self.syshandler.name))

				#Update dictionary if the batch is indicized
				if self.syshandler.exists(self.infofile):
					self.info[new_plane_set.cosmo_id][new_plane_set.geometry_id]["nbody"][str(new_plane_set.ic_index)]["plane_sets"][settings.directory_name] = dict()

		#Save a pickled copy of the settings for future reference
		with self.syshandler.open(os.path.join(new_plane_set.home_subdir,"settings.p"),"wb") as settingsfile:
			self.syshandler.pickledump(settings,settingsfile)

		#Append the name of the plane batch to a summary file
		if not self.syshandler.exists(self.infofile):
			with self.syshandler.open(os.path.join(self.home_subdir,"sets.txt"),"a") as setsfile:
				setsfile.write("{0}\n".format(settings.directory_name))

		#Return the created instance
		return new_plane_set

	####################################################################################################################################

	def getPlaneSet(self,setname):

		"""
		Instantiate a SimulationPlanes object that handles a specific plane set; returns None if the plane set does not exist

		:param setname: name of the plane set
		:type setname: str.

		:rtype: SimulationPlanes

		"""

		if not(self.syshandler.exists(os.path.join(self.storage_subdir,setname))):
			return None

		#Read plane settings from pickled file
		with self.syshandler.open(os.path.join(self.home_subdir,setname,"settings.p"),"rb") as settingsfile:
			settings = self.syshandler.pickleload(settingsfile)

		#Instantiate the SimulationPlanes object
		return SimulationPlanes(self.cosmology,self.environment,self.parameters,self.box_size,self.nside,self.ic_index,self.seed,self.ICFileBase,self.SnapshotFileBase,settings,syshandler=self.syshandler)

	####################################################################################################################################

	def linkPlaneSet(self,plane,linkname):

		"""
		Link an existing plane set to the current SimulationIC 

		:param plane: existing plane set to link
		:type plane: :py:class:`SimulationPlanes`

		:param linkname: name to give the symlinked plane set
		:type linkname: str.

		"""

		#Type check
		assert isinstance(plane,SimulationPlanes)

		#Create new plane in the directory tree
		old_plane_name = plane.settings.directory_name
		plane.settings.directory_name = linkname
		new_plane = self.newPlaneSet(plane.settings)

		#Perform the symbolic linking operation
		os.rmdir(new_plane.storage)
		os.symlink(plane.storage,new_plane.storage)
		print("[+] Linked plane set {0} to {1}".format(plane.storage,new_plane.storage))

		#Restore settings
		plane.settings.directory_name = old_plane_name

	####################################################################################################################################

	@property
	def planesets(self):

		"""
		Build a list with all the available plane sets in the current realization

		:returns: list of SimulationMaps

		"""

		plane_sets = list()
		
		try:

			if self.syshandler.exists(self.infofile):
				for set_name in self.info[self.cosmo_id][self.geometry_id]["nbody"][str(self.ic_index)]["plane_sets"].keys():
					plane_sets.append(self.getPlaneSet(set_name))
			else:
				with self.syshandler.open(os.path.join(self.home_subdir,"sets.txt"),"r") as setsfile:
					for set_name in setsfile.readlines():
						plane_sets.append(self.getPlaneSet(set_name.strip("\n")))
		
		except IOError:
			pass

		finally:
			return plane_sets


	####################################################################################################################################

	def writeNGenIC(self,settings):

		"""
		Generates the parameter file that NGenIC needs to read to generate the current initial condition

		:param settings: NGenIC tunable settings
		:type settings: NGenICSettings

		"""

		#Safety check
		assert isinstance(settings,NGenICSettings)

		#The filename is automatically generated from the class instance
		filename = os.path.join(self.home_subdir,"ngenic.param")

		#Write the parameter file
		with self.syshandler.open(filename,"w") as paramfile:

			#Mesh and grid size
			paramfile.write("Nmesh			{0}\n".format(2*self.nside))
			paramfile.write("Nsample		{0}\n".format(self.nside))

			#Box
			paramfile.write("Box 			{0:.1f}\n".format(self.box_size.to(self.kpc_over_h).value))

			#Base names for outputs
			paramfile.write("FileBase			{0}\n".format(self.ICFileBase))
			paramfile.write("OutputDir			{0}\n".format(os.path.abspath(self.ics)))

			#Glass file
			paramfile.write("GlassFile			{0}\n".format(os.path.abspath(settings.GlassFile)))

			#Tiling
			with configuration.snapshot_handler.open(os.path.abspath(settings.GlassFile)) as glass: 
				nside_glass = glass.header["num_particles_total_side"]
			paramfile.write("TileFac			{0}\n".format(self.nside//nside_glass))

			#Cosmological parameters
			paramfile.write("Omega			{0:.6f}\n".format(self.cosmology.Om0))
			paramfile.write("OmegaLambda			{0:.6f}\n".format(self.cosmology.Ode0))
			paramfile.write("OmegaBaryon			{0:.6f}\n".format(self.cosmology.Ob0))
			paramfile.write("HubbleParam			{0:.6f}\n".format(self.cosmology.h))

			if issubclass(configuration.snapshot_handler,Gadget2SnapshotDE):
				paramfile.write("w0			{0:.6f}\n".format(self.cosmology.w0))
				paramfile.write("wa			{0:.6f}\n".format(self.cosmology.wa))

			#Initial redshift
			paramfile.write("Redshift 			{0:.6f}\n".format(settings.Redshift))

			if issubclass(configuration.snapshot_handler,Gadget2SnapshotDE):
			
				#Compute the growth and velocity prefactors
				print("[+] Solving the linear growth ODE for {0}...".format(self.cosmo_id))
				g = self.cosmology.growth_factor([settings._zmaxact,settings.Redshift,0.0])

				print("[+] Computing prefactors...".format(self.cosmo_id))
				growth_prefactor = g[2,0] / g[1,0]
				vel_prefactor = -0.1 * np.sqrt(1+settings.Redshift) *(self.cosmology.H(settings.Redshift)/self.cosmology.H(0)).value * g[1,1] / g[1,0]

				paramfile.write("GrowthFactor			{0:.6f}\n".format(growth_prefactor))
				paramfile.write("VelocityPrefactor			{0:.6f}\n".format(vel_prefactor))

			#Sigma8
			paramfile.write("Sigma8				{0:.6f}\n".format(self.cosmology.sigma8))

			#Power Spectrum settings
			paramfile.write("SphereMode			{0}\n".format(settings.SphereMode))
			paramfile.write("WhichSpectrum			{0}\n".format(settings.WhichSpectrum))
			
			ngenic_ps_file = os.path.join(self.environment.home,self.cosmo_id,self.geometry_id,"ngenic_matterpower_z{0:.6f}.txt".format(0.0))

			#Check if NGen-IC power spectrum file exists, if not throw exception
			if not(self.syshandler.exists(ngenic_ps_file)) and settings.WhichSpectrum==2:
				raise IOError("NGen-IC power spectrum file {0} does not exist yet!".format(ngenic_ps_file))

			paramfile.write("FileWithInputSpectrum			{0}\n".format(ngenic_ps_file))
			
			paramfile.write("InputSpectrum_UnitLength_in_cm			{0:.6e}\n".format(settings.InputSpectrum_UnitLength_in_cm))
			paramfile.write("ReNormalizeInputSpectrum		{0}\n".format(settings.ReNormalizeInputSpectrum))
			paramfile.write("ShapeGamma			{0:.2f}\n".format(settings.ShapeGamma))
			paramfile.write("PrimordialIndex			{0:.6f}\n".format(settings.PrimordialIndex))

			#Random seed
			paramfile.write("Seed 			{0}\n".format(self.seed))

			#Files written in parallel
			paramfile.write("NumFilesWrittenInParallel			{0}\n".format(settings.NumFilesWrittenInParallel))

			#Units
			paramfile.write("UnitLength_in_cm 			{0:.6e}\n".format(settings.UnitLength_in_cm))
			paramfile.write("UnitMass_in_g			{0:.6e}\n".format(settings.UnitMass_in_g))
			paramfile.write("UnitVelocity_in_cm_per_s 			{0:.6e}\n".format(settings.UnitVelocity_in_cm_per_s))

		#Save a pickled copy of the settings for future reference
		with self.syshandler.open(os.path.join(self.home_subdir,"ngenic.p"),"wb") as settingsfile:
			self.syshandler.pickledump(settings,settingsfile)

		#Update the instance settings
		self.ngenic_settings = settings

		#Log and return
		print("[+] NGenIC parameter file {0} written on {1}".format(filename,self.syshandler.name))

	####################################################################################################################################

	def writeGadget2(self,settings):

		"""
		Generates the parameter file that Gadget2 needs to read to evolve the current initial condition in time

		:param settings: Gadget2 tunable settings
		:type settings: Gadget2Settings

		"""

		#Safety check
		assert isinstance(settings,Gadget2Settings)

		#The filename is automatically generated from the class instance
		filename = os.path.join(self.home_subdir,"gadget2.param")

		#Write the parameter file
		with self.syshandler.open(filename,"w") as paramfile:

			#File names for initial condition and outputs
			initial_condition_file = os.path.join(os.path.abspath(self.ics),self.ICFileBase)
			paramfile.write("InitCondFile			{0}\n".format(initial_condition_file))
			paramfile.write("OutputDir			{0}{1}\n".format(self.snapshots,os.path.sep))
			paramfile.write("EnergyFile			{0}\n".format(settings.EnergyFile))
			paramfile.write("InfoFile			{0}\n".format(settings.InfoFile))
			paramfile.write("TimingsFile			{0}\n".format(settings.TimingsFile))
			paramfile.write("CpuFile			{0}\n".format(settings.CpuFile))
			paramfile.write("RestartFile			{0}\n".format(settings.RestartFile))
			paramfile.write("SnapshotFileBase			{0}\n".format(self.SnapshotFileBase))

			#Use outputs in the settings to write the OutputListFilename, and set this as the output list of the code
			outputs_filename = os.path.join(self.home_subdir,"outputs.txt")
			np.savetxt(outputs_filename,settings.OutputScaleFactor)
			paramfile.write("OutputListFilename			{0}\n\n".format(os.path.abspath(outputs_filename)))

			#CPU time limit section
			paramfile.write(settings.writeSection("cpu_timings"))

			#Code options section
			paramfile.write(settings.writeSection("code_options"))

			#Initial scale factor time
			ic_filenames = self.syshandler.glob(initial_condition_file+"*")

			try:
				ic_snapshot = configuration.snapshot_handler.open(ic_filenames[0])
				paramfile.write("TimeBegin			{0}\n".format(ic_snapshot.header["scale_factor"]))
				ic_snapshot.close()
			except (IndexError,IOError):
				
				#Read the initial redshift of the simulation from the NGenIC settings
				with self.syshandler.open(os.path.join(self.home_subdir,"ngenic.p"),"rb") as ngenicfile:
					ngenic_settings = self.syshandler.pickleload(ngenicfile)
					assert isinstance(ngenic_settings,NGenICSettings)

				#Write the corresponding section of the Gadget parameter file
				paramfile.write("TimeBegin			{0:.6f}\n".format(1.0/(1+ngenic_settings.Redshift)))

			#Characteristics of run section
			paramfile.write(settings.writeSection("characteristics_of_run"))

			#Cosmological parameters
			paramfile.write("Omega0			{0:.6f}\n".format(self.cosmology.Om0))
			paramfile.write("OmegaLambda			{0:.6f}\n".format(self.cosmology.Ode0))
			paramfile.write("OmegaBaryon			{0:.6f}\n".format(self.cosmology.Ob0))
			paramfile.write("HubbleParam			{0:.6f}\n".format(self.cosmology.h))
			paramfile.write("BoxSize			{0:.6f}\n".format(self.box_size.to(self.kpc_over_h).value))

			if issubclass(configuration.snapshot_handler,Gadget2SnapshotDE):
				paramfile.write("w0			{0:.6f}\n".format(self.cosmology.w0))
				paramfile.write("wa			{0:.6f}\n\n".format(self.cosmology.wa))

			#Output frequency section
			paramfile.write(settings.writeSection("output_frequency"))

			#Accuracy of time integration section
			paramfile.write(settings.writeSection("accuracy_time_integration"))

			#Tree algorithm section
			paramfile.write(settings.writeSection("tree_algorithm"))

			#SPH section
			paramfile.write(settings.writeSection("sph"))

			#Memory allocation section
			paramfile.write(settings.writeSection("memory_allocation"))

			#System of units section
			paramfile.write(settings.writeSection("system_of_units"))

			#Softening lengths section
			paramfile.write(settings.writeSection("softening"))

		#Save a pickled copy of the settings for future reference
		with self.syshandler.open(os.path.join(self.home_subdir,"gadget2.p"),"wb") as settingsfile:
			self.syshandler.pickledump(settings,settingsfile)

		#Update the instance settings
		self.gadget_settings = settings

		#Log and exit
		print("[+] Gadget2 parameter file {0} written on {1}".format(filename,self.syshandler.name))


##############################################################
##############SimulationPlanes class##########################
##############################################################

class SimulationPlanes(SimulationIC):

	"""
	Class handler of a set of lens planes belonging to a particular simulation

	"""

	def __init__(self,cosmology,environment,parameters,box_size,nside,ic_index,seed,ICFileBase,SnapshotFileBase,settings,**kwargs):

		#Safety check
		assert isinstance(settings,PlaneSettings) or isinstance(settings,PlaneLightConeSettings)

		#Call parent constructor
		super(SimulationPlanes,self).__init__(cosmology,environment,parameters,box_size,nside,ic_index,seed,ICFileBase,SnapshotFileBase,**kwargs)
		self.settings = settings
		self.name = settings.directory_name

		#Build the name of the dedicated plane directory
		self.home_subdir = os.path.join(self.home_subdir,settings.directory_name)
		self.storage_subdir = os.path.join(self.storage_subdir,settings.directory_name)

	def __repr__(self):

		#Old representation string
		parent_repr = super(SimulationPlanes,self).__repr__().split("|")
		for i in range(2):
			parent_repr.pop(-1)

		#Count the number of plane files on disk
		planes_on_disk = self.syshandler.glob(os.path.join(self.storage_subdir,"*"))

		#Build the new representation string
		parent_repr.append("Plane set: {0} , Plane files on disk: {1}".format(self.settings.directory_name,len(planes_on_disk)))
		return " | ".join(parent_repr)

	def mkinfo(self):

		"""
		Write the plane info file in the storage subdirectory

		"""

		s = StringIO()
		info_filename = self.path("info.txt")

		#Look at all the planes present in the storage directory
		plane_files = self.ls(glob="snap*_potentialPlane0_normal0.fits")
		plane_files.sort(key=lambda pf:int(pf.split("_")[0].strip("snap")))

		#Write a line for each file
		with self.syshandler.open(info_filename,"w") as fp:
			for n,pf in enumerate(plane_files):

				#Header and snapsnot number
				header = PotentialPlane.readHeader(self.path(pf))
				nsnap = pf.split("_")[0].strip("snap")
				buf = "s={0},d={1} Mpc/h,z={2}\n".format(nsnap,header["CHI"],header["Z"])
				fp.write(buf)
				s.write(buf)

		#Report to user
		print("[+] Wrote {0} plane information lines to {1}".format(n+1,info_filename))
		s.seek(0)

		return s.read()

	###########################################################################

	#Can't call these methods anymore
	def newPlaneSet(self,*args):
		raise NotImplementedError

	@property
	def planesets(self):
		raise NotImplementedError

	def getPlaneSet(self,*args):
		raise NotImplementedError

	def writeGadget2(self,*args):
		raise NotImplementedError

	def writeNGenIC(self,*args):
		raise NotImplementedError

	def pipe_snapshots(self,*args):
		raise NotImplementedError
	
	def snapshotPath(self,*args):
		raise NotImplementedError


##############################################################
##############SimulationMaps class############################
##############################################################

class SimulationMaps(SimulationCollection):

	"""
	Class handler of a set of lensing maps, which are the final products of the lensing pipeline

	"""

	def __init__(self,cosmology,environment,parameters,box_size,nside,settings,**kwargs):

		#Safety check
		assert isinstance(settings,MapSettings)

		#Call parent constructor
		super(SimulationMaps,self).__init__(cosmology,environment,parameters,box_size,nside,**kwargs)
		self.settings = settings
		self.name = settings.directory_name

		#Build the name of the dedicated map directory
		self.home_subdir = os.path.join(self.home_subdir,settings.directory_name)
		self.storage_subdir = os.path.join(self.storage_subdir,settings.directory_name)

	def __repr__(self):

		#Count the number of map files on disk
		maps_on_disk = self.syshandler.glob(os.path.join(self.storage_subdir,"WL*"))

		#Build the new representation string
		return super(SimulationMaps,self).__repr__() + " | Map set: {0} | Map files on disk: {1} ".format(self.settings.directory_name,len(maps_on_disk))

	###########################################################################

	#Can't call these methods anymore
	def newRealization(self,seed):
		raise NotImplementedError

	@property
	def realizations(self):
		raise NotImplementedError

	def getRealization(self,ic):
		raise NotImplementedError

	def newMapSet(self,mp):
		raise NotImplementedError

	@property
	def mapsets(self):
		raise NotImplementedError

	def getMapSet(self,mp):
		raise NotImplementedError

	def newCatalog(self,ct):
		raise NotImplementedError

	@property
	def catalogs(self):
		raise NotImplementedError

	def getCatalog(self,ct):
		raise NotImplementedError

	def camb2ngenic(self,*args):
		raise NotImplementedError

	def writeCAMB(self,*args):
		raise NotImplementedError

########################################################################
##############SimulationTelescopicMaps class############################
########################################################################

class SimulationTelescopicMaps(SimulationModel):

	"""
	Class handler of a set of lensing telescopic maps, which are the final products of the lensing pipeline

	"""


	def __init__(self,cosmology,environment,parameters,collections,redshifts,settings,**kwargs):

		assert redshifts[0]==0.0,"First redshift should be 0!"
		assert len(redshifts)==len(collections)
		assert isinstance(settings,TelescopicMapSettings)

		super(SimulationTelescopicMaps,self).__init__(cosmology,environment,parameters,**kwargs)
		self.mapcollections = collections
		self.redshifts = redshifts
		self.settings = settings
		self.name = settings.directory_name

		#Build the directory names
		self.home_subdir = os.path.join(self.environment.home,self.cosmo_id,self.name)
		self.storage_subdir = os.path.join(self.environment.storage,self.cosmo_id,self.name)

	def __repr__(self):

		#Count the number of map files on disk
		maps_on_disk = self.syshandler.glob(os.path.join(self.storage_subdir,"WL*"))

		#Build the new representation string
		return super(SimulationTelescopicMaps,self).__repr__() + " | " + "-".join([c.geometry_id for c in self.mapcollections]) + " | " + "-".join([ "{0:.3f}".format(z) for z in self.redshifts ]) +" | Map files on disk: {0}".format(len(maps_on_disk))

	################################################################################################################################

	#Can't call these methods anymore
	def newCollection(self):
		raise NotImplementedError

	def getCollection(self):
		raise NotImplementedError

	@property
	def collections(self):
		raise NotImplementedError

	def newTelescopicMapSet(self):
		raise NotImplementedError

	def getTelescopicMapSet(self):
		raise NotImplementedError

	@property
	def telescopicmapsets(self):
		raise NotImplementedError		


#################################################################
##############SimulationCatalog class############################
#################################################################

class SimulationCatalog(SimulationCollection):

	"""
	Class handler of a simulated lensing catalog, which is the final products of the lensing pipeline

	"""

	def __init__(self,cosmology,environment,parameters,box_size,nside,settings,**kwargs):

		#Safety check
		assert isinstance(settings,CatalogSettings)

		#Call parent constructor
		super(SimulationCatalog,self).__init__(cosmology,environment,parameters,box_size,nside,**kwargs)
		self.settings = settings
		self.name = settings.directory_name

		#Build the name of the dedicated map directory
		self.home_subdir = os.path.join(self.home_subdir,settings.directory_name)
		self.storage_subdir = os.path.join(self.storage_subdir,settings.directory_name)

	def __repr__(self):

		#Count the number of map files on disk
		catalogs_on_disk = self.syshandler.glob(os.path.join(self.storage_subdir,"WL*"))

		#Build the new representation string
		return super(SimulationCatalog,self).__repr__() + " | Catalog set: {0} | Catalog files on disk: {1} ".format(self.settings.directory_name,len(catalogs_on_disk))

	@property 
	def subcatalogs(self):

		"""
		List the subcatalogs present in the catalog

		"""

		#List the storage subdirectory, and find the sub-catalog directories
		sub_catalog_directories = self.syshandler.glob(os.path.join(self.storage_subdir,"*-*"))
		sub_catalogs = list()

		#Build a SimulationSubCatalog instance for each sub_catalog found
		for d in sub_catalog_directories:
			basename = os.path.basename(d)
			try:
				
				#Build
				first_realization,last_realization = basename.split("-")
				first_realization = int(first_realization)
				last_realization = int(last_realization)
				sub_catalog = SimulationSubCatalog(self.cosmology,self.environment,self.parameters,self.box_size,self.nside,self.settings,syshandler=self.syshandler)
				sub_catalog.storage_subdir = os.path.join(sub_catalog.storage_subdir,basename)
				sub_catalog._first_realization = first_realization
				sub_catalog._last_realization = last_realization

				#Append
				sub_catalogs.append(sub_catalog)
			
			except ValueError:
				pass

		#Sort according to first realization and return to user
		sub_catalogs.sort(key=lambda c:c.first_realization)
		return sub_catalogs

	###########################################################################

	#Can't call these methods anymore
	def newRealization(self,seed):
		raise NotImplementedError

	@property
	def realizations(self):
		raise NotImplementedError

	def getRealization(self,ic):
		raise NotImplementedError

	def newMapSet(self,mp):
		raise NotImplementedError

	@property
	def mapsets(self):
		raise NotImplementedError

	def getMapSet(self,mp):
		raise NotImplementedError

	def newCatalog(self,ct):
		raise NotImplementedError

	@property
	def catalogs(self):
		raise NotImplementedError

	def getCatalog(self,ct):
		raise NotImplementedError

	def camb2ngenic(self,*args):
		raise NotImplementedError

	def writeCAMB(self,*args):
		raise NotImplementedError	


####################################################################################################################################


class SimulationSubCatalog(SimulationCatalog):

	"""
	Class handler of a simulated lensing sub-catalog, that contains a subset of the realizations of a bigger catalog

	"""

	def __repr__(self):

		#Count the number of map files on disk
		catalogs_on_disk = self.syshandler.glob(os.path.join(self.storage_subdir,"WL*"))

		#Build the new representation string
		return super(SimulationCatalog,self).__repr__() + " | Sub-catalog set: {0}({1}-{2}) | Catalog files on disk: {3} ".format(self.settings.directory_name,self.first_realization,self.last_realization,len(catalogs_on_disk))


	@property 
	def first_realization(self):
		return self._first_realization

	@property
	def last_realization(self):
		return self._last_realization

	@property
	def subcatalogs(self):
		raise NotImplementedError










