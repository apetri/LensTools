from __future__ import division,with_statement

import os
import re
import tarfile
import json

import numpy as np
import astropy.units as u
from astropy.cosmology import z_at_value

import configuration
from .configuration import LensToolsCosmology

from .remote import SystemHandler,LocalGit
from .settings import *

from .deploy import JobHandler
from ..simulations import Gadget2Settings,Gadget2SnapshotDE

from ..simulations.camb import CAMBSettings

#####################################################
############Parse cosmology from string##############
#####################################################

def string2cosmo(s):

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
				parameters_dict[configuration.name2attr[par]] = float(val)
		
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

	def __init__(self,batch):
		
		self.batch = batch
		self.dictionary = batch.info

	def __enter__(self):
		return self

	def __exit__(self,type,value,tb):
		pass

	def update(self):
		dictionary_file = os.path.join(self.batch.home_subdir,configuration.json_tree_file)
		with self.batch.syshandler.open(dictionary_file,"w") as fp:
			fp.write(json.dumps(self.dictionary))


#####################################################
##############SimulationBatch class##################
#####################################################

class SimulationBatch(object):

	"""
	Class handler of a batch of weak lensing simulations that share the same environment settings

	"""

	@classmethod
	def current(cls,syshandler=configuration.syshandler,indicize=False):

		"""
		This method looks in the current directory and looks for a configuration file named "environment.ini"; if it finds one, it returns a SimulationBatch instance that corresponds to the one pointed to by "environment.ini"

		:param syshandler: system handler that allows to override the methods used to create directories and do I/O from files, must implement the abstract type SystemHandler
		:type syshandler: SystemHandler

		:returns: Simulation batch pointed to by "environment.ini", or None
		:rtype: SimulationBatch

		"""

		if not(os.path.exists("environment.ini")):
			return None

		env = EnvironmentSettings.read("environment.ini")
		return cls(env,syshandler,indicize)

	@property 
	def home_subdir(self):
		return self.environment.home

	@property
	def storage_subdir(self):
		return self.environment.storage
		

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
			with InfoDict(self) as info:
				info.update()

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
		dirnames = [ os.path.basename(n) for n in self.syshandler.glob(os.path.join(self.environment.home,"*")) ]

		#Decide which of the directories actually correspond to cosmological models
		for dirname in dirnames:
			cosmo_parsed = string2cosmo(dirname)
			if cosmo_parsed is not None:
				models.append(SimulationModel(cosmology=cosmo_parsed[0],environment=self.environment,parameters=cosmo_parsed[1],syshandler=self.syshandler))

		#Return the list with the available models
		return models

	#Alias for available models
	@property
	def models(self):
		return self.available

	##############################################################################################################################

	@property
	def info(self):

		"""
		Returns summary info of the simulation batch corresponding to the current environment

		:returns: info in dictionary format

		"""

		#Load info from tree file if available
		tree_file_path = os.path.join(self.home_subdir,configuration.json_tree_file)
		if self.syshandler.exists(tree_file_path):
			with self.syshandler.open(tree_file_path,"r") as fp:
				info_dict = json.loads(fp.read())
			return info_dict

		#Information will be returned in dictionary format
		info_dict = dict()

		#Start with the available models
		available_models = self.available
		for model in available_models:
			info_dict[model.cosmo_id] = dict()

			#Follow with the collections 
			for collection in model.collections:
				info_dict[model.cosmo_id][collection.geometry_id] = dict()
				info_dict[model.cosmo_id][collection.geometry_id]["map_sets"] = dict()
				info_dict[model.cosmo_id][collection.geometry_id]["catalogs"] = dict()

				#Check if there are any map sets or catalogs present
				try:
					
					with self.syshandler.open(os.path.join(collection.home_subdir,"sets.txt"),"r") as setsfile:
						
						for line in setsfile.readlines():
							if line=="":
								continue
							map_set = collection.getMapSet(line.strip("\n"))
							maps_on_disk = self.syshandler.glob(os.path.join(map_set.storage_subdir,"WL*"))
							info_dict[model.cosmo_id][collection.geometry_id]["map_sets"][map_set.settings.directory_name] = dict() 
							info_dict[model.cosmo_id][collection.geometry_id]["map_sets"][map_set.settings.directory_name]["num_maps"] = len(maps_on_disk)

				except IOError:
					pass

				
				try:
					
					with self.syshandler.open(os.path.join(collection.home_subdir,"catalogs.txt"),"r") as setsfile:
						
						for line in setsfile.readlines():
							if line=="":
								continue
							
							catalog = collection.getCatalog(line.strip("\n"))
							catalogs_on_disk = self.syshandler.glob(os.path.join(catalog.storage_subdir,"WL*"))
							info_dict[model.cosmo_id][collection.geometry_id]["catalogs"][catalog.settings.directory_name] = dict() 
							info_dict[model.cosmo_id][collection.geometry_id]["catalogs"][catalog.settings.directory_name]["num_catalogs"] = len(catalogs_on_disk)

				except IOError:
					pass

				#Follow with the realizations
				for r in collection.realizations:
					
					info_dict[model.cosmo_id][collection.geometry_id][r.ic_index] = dict()
					info_dict[model.cosmo_id][collection.geometry_id][r.ic_index]["plane_sets"] = dict()

					#Make number of ics and snapshots on disk available to user
					ics_on_disk = self.syshandler.glob(os.path.join(r.ics_subdir,r.ICFilebase+"*"))
					snap_on_disk = self.syshandler.glob(os.path.join(r.snapshot_subdir,r.SnapshotFileBase+"*"))
					info_dict[model.cosmo_id][collection.geometry_id][r.ic_index]["ics"] =  len(ics_on_disk)
					info_dict[model.cosmo_id][collection.geometry_id][r.ic_index]["snapshots"] = len(snap_on_disk)

					#Check if there are any plane sets present
					try:

						with self.syshandler.open(os.path.join(r.home_subdir,"sets.txt"),"r") as setsfile:

							for line in setsfile.readlines():
								if line=="":
									continue
								
								plane_set = r.getPlaneSet(line.strip("\n"))
								planes_on_disk = self.syshandler.glob(os.path.join(plane_set.storage_subdir,"*Plane*"))
								info_dict[model.cosmo_id][collection.geometry_id][r.ic_index]["plane_sets"][plane_set.settings.directory_name] = dict()
								info_dict[model.cosmo_id][collection.geometry_id][r.ic_index]["plane_sets"][plane_set.settings.directory_name]["num_planes"] = len(planes_on_disk)

					except IOError:
						pass


		#Return to user
		return info_dict


	##############################################################################################################################################

	def list(self,resource=None,which=None,chunk_size=10,**kwargs):

		"""
		Lists the available resources in the simulation batch (collections,mapsets,etc...)

		:param resource: custom function to call on each batch.available element, must return a string. If None the list of Storage model directories is returned
		:type resource: None or callable

		:param which: extremes of the model numbers to get (if None all models are processed)
		:type which: tuple. 

		:param chunk_size: size of output chunk
		:type chunk_size: int.

		:param kwargs: the keyword arguments are passed to resource
		:type kwargs: dict.

		:returns: requested resources
		:rtype: list.

		"""

		#Available models
		if which is None:
			models = self.available
		else:
			models = self.available.__getslice__(*which)

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

	def archive(self,name,pool=None,**kwargs):

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
		resource_chunks = self.list(**kwargs)

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
		cosmo_parsed = string2cosmo(cosmo_id)

		#Return the SimulationModel instance
		if cosmo_parsed is not None:
			return SimulationModel(cosmology=cosmo_parsed[0],environment=self.environment,parameters=cosmo_parsed[1],syshandler=self.syshandler)
		else:
			return None

	###########################################################################################################################################
	##########################################Job submission scripts###########################################################################
	###########################################################################################################################################

	def writeCAMBSubmission(self,realization_list,job_settings,job_handler,chunks=1,**kwargs):

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

				parameter_file = os.path.join(collection.home_subdir,"camb.param")
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

	def writeNGenICSubmission(self,realization_list,job_settings,job_handler,chunks=1,**kwargs):

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
		
			#Arguments for the executable
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

				parameter_file = os.path.join(r.home_subdir,"ngenic.param")
				if not(self.syshandler.exists(parameter_file)):
					raise IOError("NGenIC parameter file at {0} does not exist yet!".format(parameter_file))

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

			#Override settings
			job_settings.num_cores = job_settings.cores_per_simulation

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


	def writeGadget2Submission(self,realization_list,job_settings,job_handler,chunks=1,**kwargs):
		
		"""
		Writes Gadget2 submission script

		:param realization_list: list of ics to generate in the form "cosmo_id|geometry_id|icN"
		:type realization_list: list. of str.

		:param chunks: number of independent jobs in which to split the submission (one script per job will be written)
		:type chunks: int. 

		:param job_settings: settings for the job (resources, etc...)
		:type job_settings: JobSettings

		:param job_handler: handler of the cluster specific features (job scheduler, architecture, etc...)
		:type job_handler: JobHandler

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

				parameter_file = os.path.join(r.home_subdir,"gadget2.param")
				if not(self.syshandler.exists(parameter_file)):
					raise IOError("Gadget2 parameter file at {0} does not exist yet!".format(parameter_file))

				executables.append(job_settings.path_to_executable + " " + "{0} {1} {2}".format(1,job_settings.cores_per_simulation,parameter_file))

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
				with self.syshandler.open(os.path.join(r.home_subdir,"gadget2.p"),"r") as settingsfile:
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
						assert raytracing_settings.lens_map_realizations%job_settings.cores_per_simulation==0,"The number of map realizations must be a multiple of the number of cores per simulation!"
					except AssertionError:
						raytracing_settings = CatalogSettings.read(config_file)
						assert raytracing_settings.lens_catalog_realizations%job_settings.cores_per_simulation==0,"The number of map realizations must be a multiple of the number of cores per simulation!"

				elif len(parts)==1:
					raytracing_settings = TelescopicMapSettings.read(config_file)
					assert raytracing_settings.lens_map_realizations%job_settings.cores_per_simulation==0,"The number of map realizations must be a multiple of the number of cores per simulation!"
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




#####################################################
##############SimulationModel class##################
#####################################################

class SimulationModel(object):

	"""
	Class handler of a weak lensing simulation model, defined by a set of cosmological parameters

	"""

	@property
	def cosmo_id(self):
		base = "{0}{1:."+str(configuration.cosmo_id_digits)+"f}"
		return "_".join([ base.format(p,getattr(self.cosmology,configuration.name2attr[p])) for p in self.parameters if (hasattr(self.cosmology,configuration.name2attr[p]) and getattr(self.cosmology,configuration.name2attr[p]) is not None)])

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
			representation_parameters.append(("{0}={1:."+str(configuration.cosmo_id_digits)+"f}").format(p,getattr(self.cosmology,configuration.name2attr[p])))

		return "<"+ " , ".join(representation_parameters) + ">"


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

		full_path = os.path.join(getattr(self,where),filename)
		if not(self.syshandler.exists(full_path)):
			return None

		return full_path

	################################################################################################################################

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
		return [os.path.basename(f) for f in self.syshandler.glob(os.path.join(search_path,glob))] 

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
		if type(box_size)==str:
			parts = box_size.split("b")
			nside = int(parts[0])
			box_size = float(parts[1]) * self.Mpc_over_h

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

		collection_names = [ os.path.basename(d) for d in self.syshandler.glob(os.path.join(self.home_subdir,"*")) ]
		collection_list = list()

		for name in collection_names:
			
			try:
				nside,box = name.split("b")
				collection_list.append(SimulationCollection(self.cosmology,self.environment,self.parameters,box_size=float(box)*self.Mpc_over_h,nside=int(nside),syshandler=self.syshandler))
			except ValueError:
				pass

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
		with self.syshandler.open(os.path.join(map_set.home_subdir,"settings.p"),"w") as settingsfile:
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
				with self.syshandler.open(os.path.join(self.home_subdir,name,"settings.p"),"r") as settingsfile:
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
			with self.syshandler.open(os.path.join(self.home_subdir,name,"settings.p"),"r") as settingsfile:
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


	def newCollection(self):
		raise TypeError("This method should be called on SimulationModel instances!")

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

	def newRealization(self,seed=0):
		
		#Check if there are already generated initial conditions in there
		ics_present = self.syshandler.glob(os.path.join(self.storage_subdir,"ic*"))
		new_ic_index = len(ics_present) + 1

		#Generate the new initial condition
		newIC = SimulationIC(self.cosmology,self.environment,self.parameters,self.box_size,self.nside,new_ic_index,seed,syshandler=self.syshandler)

		#Make dedicated directories for new initial condition,ics and snapshots
		for d in [newIC.home_subdir,newIC.storage_subdir,newIC.ics_subdir,newIC.snapshot_subdir]:

			if not self.syshandler.exists(d):
				self.syshandler.mkdir(d)
				print("[+] {0} created on {1}".format(d,self.syshandler.name))

		#Make new file with the number of the seed
		with self.syshandler.open(os.path.join(newIC.home_subdir,"seed"+str(seed)),"w") as seedfile:
			pass

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
		seed_filename = os.path.basename(self.syshandler.glob(os.path.join(newIC.home_subdir,"seed*"))[0])
		newIC.seed = int(seed_filename.strip("seed"))

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

		ic_list = list()
		ic_numbers = [ os.path.basename(d).strip("ic") for d in self.syshandler.glob(os.path.join(self.home_subdir,"ic*")) ]

		for ic in ic_numbers:

			seed = int(os.path.basename(self.syshandler.glob(os.path.join(self.home_subdir,"ic"+ic,"seed*"))[0]).strip("seed"))
			ic_list.append(SimulationIC(self.cosmology,self.environment,self.parameters,self.box_size,self.nside,int(ic),seed,syshandler=self.syshandler))

		#Sort according to ic_index
		ic_list.sort(key=lambda r:r.ic_index)

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

		#Save a picked copy of the settings to use for future reference
		with self.syshandler.open(os.path.join(map_set.home_subdir,"settings.p"),"w") as settingsfile:
			self.syshandler.pickledump(settings,settingsfile)

		#Append the name of the map batch to a summary file
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
		with self.syshandler.open(os.path.join(self.home_subdir,setname,"settings.p"),"r") as settingsfile:
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

		#Save a picked copy of the settings to use for future reference
		with self.syshandler.open(os.path.join(catalog.home_subdir,"settings.p"),"w") as settingsfile:
			self.syshandler.pickledump(settings,settingsfile)

		#Append the name of the map batch to a summary file
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
		with self.syshandler.open(os.path.join(self.home_subdir,catalog_name,"settings.p"),"r") as settingsfile:
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
			with self.syshandler.open(os.path.join(self.home_subdir,"catalogs.txt"),"r") as setsfile:
				for catalog_name in setsfile.readlines():
					catalogs.append(self.getCatalog(catalog_name.strip("\n")))
		except IOError:
			pass
		finally:
			return catalogs


	#################################################################################################################################

	def writeCAMB(self,z,settings):

		"""
		Generates the parameter file that CAMB needs to read to evolve the current initial condition in time

		:param settings: CAMB tunable settings
		:type settings: CAMBSettings

		:param z: redshifts at which CAMB needs to compute the matter power spectrum
		:type z: array.

		"""

		#Safety type check
		assert isinstance(settings,CAMBSettings)
		if type(z)==np.float:
			z = np.array([z])

		#Write the parameter file
		camb_filename = os.path.join(self.home_subdir,"camb.param") 
		with self.syshandler.open(camb_filename,"w") as paramfile:
			paramfile.write(settings.write(output_root=os.path.join(self.home_subdir,"camb"),cosmology=self.cosmology,redshifts=z))

		print("[+] {0} written on {1}".format(camb_filename,self.syshandler.name))

		#Save a pickled copy of the settings
		with self.syshandler.open(os.path.join(self.home_subdir,"camb.p"),"w") as settingsfile:
			self.syshandler.pickledump(settings,settingsfile)


	################################################################################################################################

	def camb2ngenic(self,z):

		"""
		Read CAMB power spectrum file and convert it in a N-GenIC readable format

		:param z: redshift of the matter power spectrum file to convert
		:type z: float.

		"""

		camb_ps_file = os.path.join(self.environment.home,self.cosmo_id,self.geometry_id,"camb_matterpower_z{0:.6f}.dat".format(z))
		if not(self.syshandler.exists(camb_ps_file)):
			raise IOError("CAMB matter power spectrum file {0} does not exist yet!!".format(camb_ps_file))

		k,Pk = np.loadtxt(camb_ps_file,unpack=True)
		lgk,lgP = _camb2ngenic(k,Pk)

		ngenic_ps_file = os.path.join(self.environment.home,self.cosmo_id,self.geometry_id,"ngenic_matterpower_z{0:.6f}.txt".format(z))
		np.savetxt(ngenic_ps_file,np.array([lgk,lgP]).T)

		print("[+] CAMB matter power spectrum at {0} converted into N-GenIC readable format at {1}".format(camb_ps_file,ngenic_ps_file))



##########################################################
##############SimulationIC class##########################
##########################################################

class SimulationIC(SimulationCollection):

	"""
	Class handler of a simulation with a defined initial condition

	"""

	def __init__(self,cosmology,environment,parameters,box_size,nside,ic_index,seed,ICFilebase="ics",SnapshotFileBase="snapshot",**kwargs):

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
		self.ICFilebase = ICFilebase
		self.SnapshotFileBase = SnapshotFileBase

		#Try to load in the simulation settings, if any are present
		try:
			with self.syshandler.open(os.path.join(self.home_subdir,"ngenic.p"),"r") as settingsfile:
				self.ngenic_settings = self.syshandler.pickleload(settingsfile)
		except IOError:
			pass

		try:
			with self.syshandler.open(os.path.join(self.home_subdir,"gadget2.p"),"r") as settingsfile:
				self.gadget_settings = self.syshandler.pickleload(settingsfile)
		except IOError:
			pass

	def __repr__(self):

		#Check if snapshots and/or initial conditions are present
		ics_on_disk = self.syshandler.glob(os.path.join(self.ics_subdir,self.ICFilebase+"*"))
		snap_on_disk = self.syshandler.glob(os.path.join(self.snapshot_subdir,self.SnapshotFileBase+"*"))

		return super(SimulationIC,self).__repr__() + " | ic={0},seed={1} | IC files on disk: {2} | Snapshot files on disk: {3}".format(self.ic_index,self.seed,len(ics_on_disk),len(snap_on_disk))

	def newRealization(self,seed):
		raise TypeError("This method should be called on SimulationCollection instances!")

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

	def newPlaneSet(self,settings):

		#Safety check
		assert isinstance(settings,PlaneSettings)

		#Instantiate a SimulationPlanes object
		new_plane_set = SimulationPlanes(self.cosmology,self.environment,self.parameters,self.box_size,self.nside,self.ic_index,self.seed,self.ICFilebase,self.SnapshotFileBase,settings,syshandler=self.syshandler)

		#Create the dedicated directory if not present already
		for d in [new_plane_set.home_subdir,new_plane_set.storage_subdir]:
			if not(self.syshandler.exists(d)):
				self.syshandler.mkdir(d)
				print("[+] {0} created on {1}".format(d,self.syshandler.name))

		#Save a pickled copy of the settings for future reference
		with self.syshandler.open(os.path.join(new_plane_set.home_subdir,"settings.p"),"w") as settingsfile:
			self.syshandler.pickledump(settings,settingsfile)

		#Append the name of the plane batch to a summary file
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
		with self.syshandler.open(os.path.join(self.home_subdir,setname,"settings.p"),"r") as settingsfile:
			settings = self.syshandler.pickleload(settingsfile)

		#Instantiate the SimulationPlanes object
		return SimulationPlanes(self.cosmology,self.environment,self.parameters,self.box_size,self.nside,self.ic_index,self.seed,self.ICFilebase,self.SnapshotFileBase,settings,syshandler=self.syshandler)

	####################################################################################################################################

	@property
	def planesets(self):

		"""
		Build a list with all the available plane sets in the current realization

		:returns: list of SimulationMaps

		"""

		plane_sets = list()
		
		try:
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
			paramfile.write("FileBase			{0}\n".format(self.ICFilebase))
			paramfile.write("OutputDir			{0}\n".format(os.path.abspath(self.ics_subdir)))

			#Glass file
			paramfile.write("GlassFile			{0}\n".format(os.path.abspath(settings.GlassFile)))

			#Tiling
			glass = Gadget2SnapshotDE.open(os.path.abspath(settings.GlassFile))
			nside_glass = glass.header["num_particles_total_side"]
			glass.close()
			paramfile.write("TileFac			{0}\n".format(self.nside//nside_glass))

			#Cosmological parameters
			paramfile.write("Omega			{0:.6f}\n".format(self.cosmology.Om0))
			paramfile.write("OmegaLambda			{0:.6f}\n".format(self.cosmology.Ode0))
			paramfile.write("OmegaBaryon			{0:.6f}\n".format(self.cosmology.Ob0))
			paramfile.write("HubbleParam			{0:.6f}\n".format(self.cosmology.h))
			paramfile.write("w0			{0:.6f}\n".format(self.cosmology.w0))
			paramfile.write("wa			{0:.6f}\n".format(self.cosmology.wa))

			#Initial redshift
			paramfile.write("Redshift 			{0:.6f}\n".format(settings.Redshift))

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
		with self.syshandler.open(os.path.join(self.home_subdir,"ngenic.p"),"w") as settingsfile:
			self.syshandler.pickledump(settings,settingsfile)

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
			initial_condition_file = os.path.join(os.path.abspath(self.ics_subdir),self.ICFilebase)
			paramfile.write("InitCondFile			{0}\n".format(initial_condition_file))
			paramfile.write("OutputDir			{0}{1}\n".format(self.snapshot_subdir,os.path.sep))
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
				ic_snapshot = Gadget2SnapshotDE.open(ic_filenames[0])
				paramfile.write("TimeBegin			{0}\n".format(ic_snapshot.header["scale_factor"]))
				ic_snapshot.close()
			except (IndexError,IOError):
				
				#Read the initial redshift of the simulation from the NGenIC settings
				with self.syshandler.open(os.path.join(self.home_subdir,"ngenic.p"),"r") as ngenicfile:
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
		with self.syshandler.open(os.path.join(self.home_subdir,"gadget2.p"),"w") as settingsfile:
			self.syshandler.pickledump(settings,settingsfile)

		#Log and exit
		print("[+] Gadget2 parameter file {0} written on {1}".format(filename,self.syshandler.name))


##############################################################
##############SimulationPlanes class##########################
##############################################################

class SimulationPlanes(SimulationIC):

	"""
	Class handler of a set of lens planes belonging to a particular simulation

	"""

	def __init__(self,cosmology,environment,parameters,box_size,nside,ic_index,seed,ICFilebase,SnapshotFileBase,settings,**kwargs):

		#Safety check
		assert isinstance(settings,PlaneSettings)

		#Call parent constructor
		super(SimulationPlanes,self).__init__(cosmology,environment,parameters,box_size,nside,ic_index,seed,ICFilebase,SnapshotFileBase,**kwargs)
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










