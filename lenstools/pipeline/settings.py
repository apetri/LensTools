import os
import ast
import StringIO

from ConfigParser import NoOptionError
import json

import numpy as np
import astropy.units as u

from ..simulations.settings import select_parser,LTSettings
from ..simulations.camb import CAMBSettings


############################################################
#############EnvironmentSettings class######################
############################################################

class EnvironmentSettings(LTSettings):

	"""
	This class handles the system specific environment settings, such as directory paths, modules, etc...

	"""

	def __init__(self,home="SimTest/Home",storage="SimTest/Storage"):

		"""
		Creates the home (meant to store small files like execution scripts) and storage (meant to store large files like simulation outputs) directories

		:param home: name of the simulation home directory
		:type home: str.

		:param storage: name of the simulation mass storage directory
		:type storage: str.

		"""

		self.home = home
		self.storage = storage
		self.cosmo_id_digits = 3
		self.name2attr = {"Om":"Om0","Ol":"Ode0","w":"w0","wa":"wa","h":"h","Ob":"Ob0","si":"sigma8","ns":"ns"}
		self.json_tree_file = ".tree.json"


	@classmethod
	def get(cls,options):

		#Check that the config file has the appropriate section
		section = "EnvironmentSettings"
		assert options.has_section(section),"No {0} section in configuration file {1}".format(section,options.filename)

		#Fill in the appropriate fields and return to user
		settings = cls(home=options.get(section,"home"),storage=options.get(section,"storage"))
		
		try:
			name2attr = options.get(section,"name2attr")
			if type(name2attr)==dict:
				settings.name2attr = name2attr
			else:
				settings.name2attr = json.loads(name2attr)
		except NoOptionError:
			pass

		try:
			settings.cosmo_id_digits = options.getint(section,"cosmo_id_digits")
		except NoOptionError:
			pass

		try:
			settings.json_tree_file = options.get(section,"json_tree_file")
		except NoOptionError:
			pass

		#Return
		return settings

#################################################
###########NGenICSettings class##################
#################################################

class NGenICSettings(LTSettings):

	"""
	Class handler of NGenIC settings
	
	"""

	def __init__(self,**kwargs):

		self.GlassFile = "dummy_glass_little_endian.dat"
		self.Redshift = 100.0

		self.SphereMode = 1 
		self.WhichSpectrum = 2

		self.InputSpectrum_UnitLength_in_cm = 3.085678e24
		self.ReNormalizeInputSpectrum = 1
		self.ShapeGamma = 0.21
		self.PrimordialIndex = 1.0

		self.NumFilesWrittenInParallel = 4

		self.UnitLength_in_cm = 3.085678e21
		self.UnitMass_in_g = 1.989e43
		self.UnitVelocity_in_cm_per_s = 1e5

		#Typically do not touch these, needed for prefactors calculations
		self._zmaxact = 1000.0
		self._iwmode = 3

		#Allow for kwargs override
		for key in kwargs.keys():
			setattr(self,key,kwargs[key])

############################################################
################Gadget2Settings class#######################
############################################################

class Gadget2Settings(LTSettings):

	"""
	Class handler of the tunable settings in a Gadget2 run

	"""

	file_names = ["InitCondFile","OutputDir","EnergyFile","InfoFile","TimingsFile","CpuFile","RestartFile","SnapshotFileBase","OutputListFilename"]
	cpu_timings = ["TimeLimitCPU","ResubmitOn","ResubmitCommand"]
	code_options = ["ICFormat","SnapFormat","ComovingIntegrationOn","TypeOfTimestepCriterion","OutputListOn","PeriodicBoundariesOn"]
	characteristics_of_run = ["TimeMax"]
	output_frequency = ["TimeBetSnapshot","TimeOfFirstSnapshot","CpuTimeBetRestartFile","TimeBetStatistics","NumFilesPerSnapshot","NumFilesWrittenInParallel"]
	accuracy_time_integration = ["ErrTolIntAccuracy","MaxRMSDisplacementFac","CourantFac","MaxSizeTimestep","MinSizeTimestep"]
	tree_algorithm = ["ErrTolTheta","TypeOfOpeningCriterion","ErrTolForceAcc","TreeDomainUpdateFrequency"]
	sph = ["DesNumNgb","MaxNumNgbDeviation","ArtBulkViscConst","InitGasTemp","MinGasTemp"]
	memory_allocation = ["PartAllocFactor","TreeAllocFactor","BufferSize"]
	system_of_units = ["UnitLength_in_cm","UnitMass_in_g","UnitVelocity_in_cm_per_s","GravityConstantInternal"]
	softening = ["MinGasHsmlFractional","SofteningGas","SofteningHalo","SofteningDisk","SofteningBulge","SofteningStars","SofteningBndry","SofteningGasMaxPhys","SofteningHaloMaxPhys","SofteningDiskMaxPhys","SofteningBulgeMaxPhys","SofteningStarsMaxPhys","SofteningBndryMaxPhys"]

	def __init__(self,**kwargs):

		#Default outputs
		self.OutputScaleFactor = np.array([0.2463687286034,0.25331915596808,0.2603792960781,0.26755062217624,0.27483471475984,0.282233268286,0.28974809827109,0.29738114881152,0.30513450055509,0.31301037915502,0.32101116424033,0.32913939894726,0.3373978000372,0.34578926866836,0.35431690185511,0.36298400467612,0.37179410329025,0.38075095882651,0.38985858222059,0.39912125007794,0.40854352165158,0.4181302570317,0.42788663665433,0.43781818224776,0.4479307793476,0.45823070152614,0.46872463649679,0.47941971427274,0.49032353757851,0.50144421473605,0.51279039527177,0.52437130852031,0.53619680553253,0.54827740463258,0.56062434101017,0.57324962078195,0.58616608000982,0.59938744922579,0.61292842408364,0.62680474283878,0.64103327145057,0.65563209720868,0.67062063190864,0.68601972574487,0.70185179325524,0.71814095284413,0.73491318163551,0.75219648767026,0.77002110176951,0.78841969174729,0.80742760208188,0.82708312265889,0.84742779079644,0.86850673147378,0.89036904153293,0.91306822464037,0.9366626850185,0.96121628943351,0.98679900871668,1.0])

		#File names
		self.InitCondFile = "gadget_ic"
		self.OutputDir = "snapshots"
		self.EnergyFile = "energy.txt"
		self.InfoFile = "info.txt"
		self.TimingsFile = "timings.txt"
		self.CpuFile = "cpu.txt"
		self.RestartFile = "restart"
		self.SnapshotFileBase = "snapshot"
		self.OutputListFilename = "outputs.txt"

		#CPU Timings
		self.TimeLimitCPU = 1.0*u.day
		self.ResubmitOn = 0
		self.ResubmitCommand = "my-scriptfile"

		#Code options
		self.ICFormat  = 1
		self.SnapFormat = 1
		self.ComovingIntegrationOn = 1
		self.TypeOfTimestepCriterion = 0
		self.OutputListOn = 1
		self.PeriodicBoundariesOn = 1

		#Caracteristics of run  
		self.TimeMax = 1.0

		#Output frequency
		self.TimeBetSnapshot = 0.5
		self.TimeOfFirstSnapshot = 0
		self.CpuTimeBetRestartFile = 12.5*u.hour 
		self.TimeBetStatistics = 0.05
		self.NumFilesPerSnapshot = 16
		self.NumFilesWrittenInParallel = 8

		#Accuracy of time integration
		self.ErrTolIntAccuracy = 0.025 
		self.MaxRMSDisplacementFac = 0.2
		self.CourantFac = 0.15     
		self.MaxSizeTimestep = 0.02
		self.MinSizeTimestep = 0.0


		#Tree algorithm, force accuracy, domain update frequency
		self.ErrTolTheta = 0.45
		self.TypeOfOpeningCriterion = 1
		self.ErrTolForceAcc = 0.005
		self.TreeDomainUpdateFrequency = 0.025

		
		#Further parameters of SPH
		self.DesNumNgb = 33
		self.MaxNumNgbDeviation = 2
		self.ArtBulkViscConst = 0.8
		self.InitGasTemp = 1000.0    
		self.MinGasTemp = 50.0    

		#Memory allocation
		self.PartAllocFactor = 1.3    
		self.TreeAllocFactor = 0.7
		self.BufferSize = 20*u.Mbyte 


		#System of units
		self.UnitLength_in_cm = 3.085678e21       # ;  1.0 kpc 
		self.UnitMass_in_g = 1.989e43    #;  1.0e10 solar masses 
		self.UnitVelocity_in_cm_per_s = 1.0e5  # ;  1 km/sec 
		self.GravityConstantInternal = 0


		#Softening lengths
		self.MinGasHsmlFractional = 0.25
		self.SofteningGas = 0
		self.SofteningHalo = 9.000000
		self.SofteningDisk = 0
		self.SofteningBulge = 0
		self.SofteningStars = 0
		self.SofteningBndry = 0

		self.SofteningGasMaxPhys = 0
		self.SofteningHaloMaxPhys = 9.000000
		self.SofteningDiskMaxPhys = 0
		self.SofteningBulgeMaxPhys = 0
		self.SofteningStarsMaxPhys = 0
		self.SofteningBndryMaxPhys = 0

		#Allow for kwargs override
		for key in kwargs.keys():
			setattr(self,key,kwargs[key])

	def sections(self):

		return [ "file_names","cpu_timings","code_options","characteristics_of_run","output_frequency","accuracy_time_integration","tree_algorithm","sph","memory_allocation","system_of_units","softening" ]

	def showSection(self,section):

		if section not in self.sections():
			raise ValueError("Parameter file does not admit a section named {0}".format(section))

		for option in getattr(self,section):
			print("{0} = {1}".format(option,getattr(self,option)))

	def show(self):

		for section in self.sections():
			print(section+":\n")
			self.showSection(section)
			print("\n")	

	def writeSection(self,section):

		"""
		Writes the corresponding section of the Gadget2 parameter file

		"""

		output = StringIO.StringIO()

		#Write preamble
		output.write("% {0}\n\n".format(section))

		#Cycle through options
		for option in getattr(self,section):

			#Read the corresponding value
			value = getattr(self,option)

			#Convert units as necessary
			if type(value)==u.quantity.Quantity:
				
				if value.unit.physical_type=="time":
					value = value.to(u.s).value
				elif value.unit.physical_type=="speed":
					value = value.to(u.cm/u.s).value
				elif "byte" in value.unit.to_string():
					value = value.to(u.Mbyte).value

			#Write the line
			output.write("{0}		{1}\n".format(option,value))

		#Finish
		output.write("\n\n")
		output.seek(0)

		return output.read()
		

	@classmethod
	def default(cls):

		"""
		Generate default settings
		"""

		return cls()


#################################################
###########PlaneSettings class###################
#################################################

class PlaneSettings(LTSettings):

	"""
	Class handler of plane generation settings

	"""

	def __init__(self,**kwargs):

		#Name of the planes batch
		self.directory_name = "Planes"
		
		#Use the pickled options generated at the moment of the batch generation (advised)
		self.override_with_local = True

		self.format = "fits"
		self.name_format = "snap{0}_{1}Plane{2}_normal{3}.{4}"

		self.plane_resolution = 128
		self.first_snapshot = 46
		self.last_snapshot = 58
		self.snapshots = None
		self.cut_points = np.array([7.5/0.7]) * u.Mpc
		self.thickness = (2.5/0.7) * u.Mpc 
		self.length_unit = u.Mpc
		self.normals = range(3)

		#Optional, not usually changed
		self.thickness_resolution = 1
		self.smooth = 1
		self.kind = "potential"

		#Allow for kwargs override
		for key in kwargs.keys():
			setattr(self,key,kwargs[key])

	@classmethod
	def get(cls,options):

		#Check that the config file has the appropriate section
		section = "PlaneSettings"
		assert options.has_section(section),"No {0} section in configuration file {1}".format(section,options.filename)

		#Fill in the appropriate fields
		settings = cls()
		
		settings.directory_name = options.get(section,"directory_name")
		settings.override_with_local = options.getboolean(section,"override_with_local")
		
		settings.format = options.get(section,"format")
		try:
			settings.name_format = options.get(section,"name_format")
		except NoOptionError:
			pass
		
		settings.plane_resolution = options.getint(section,"plane_resolution")

		#Snapshots to process: either bounds (first,last) or snapshot list
		try:
			settings.first_snapshot = options.getint(section,"first_snapshot")
			settings.last_snapshot = options.getint(section,"last_snapshot")
		except ValueError:
			settings.first_snapshot = None
			settings.last_snapshot = None

		try:
			snapshots = options.get(section,"snapshots") 
			settings.snapshots = [ int(n) for n in snapshots.split(",") ]
		except (NoOptionError,ValueError):
			settings.snapshots = None

		#Check that either a bound specification or a list were provided, not both
		if not(((settings.first_snapshot is not None) and (settings.last_snapshot is not None))^(settings.snapshots is not None)):
			raise ValueError("You must specify one, and only one, between (first_snapshot,last_snapshot) or a snapshot list!")

		#Length units
		settings.length_unit = getattr(u,options.get(section,"length_unit"))

		#Cut points
		settings.cut_points = np.array([ float(p) for p in options.get(section,"cut_points").split(",") ]) * settings.length_unit
		settings.thickness = options.getfloat(section,"thickness") * settings.length_unit

		#Normals
		settings.normals = [ int(n) for n in options.get(section,"normals").split(",") ]

		#Optionals
		try:
			settings.thickness_resolution = options.getint(section,"thickness_resolution")
		except NoOptionError:
			pass

		try:
			settings.smooth = options.getint(section,"smooth")
		except NoOptionError:
			pass

		try:
			settings.kind = options.get(section,"kind")
		except NoOptionError:
			pass

		#Return to user
		return settings


#################################################
###########MapSettings class#####################
#################################################

class MapSettings(LTSettings):

	"""
	Class handler of map generation settings

	"""

	_section = "MapSettings"

	def __init__(self,**kwargs):

		self._init_commmon()
		self._init_plane_set()
		self._init_randomizer()

		#Allow for kwargs override
		for key in kwargs.keys():
			setattr(self,key,kwargs[key])

	def _init_commmon(self):

		#Names of the map batch
		self.directory_name = "Maps"

		#Use the options generated at the moment of the batch generation (advised)
		self.override_with_local = True

		self.format = "fits"
		self.plane_format = "fits"
		self.plane_name_format = "snap{0}_potentialPlane{1}_normal{2}.{3}"

		self.map_resolution = 128
		self.map_angle = 1.6 * u.deg
		self.angle_unit = u.deg
		self.source_redshift = 2.0

		#Random seed used to generate multiple map realizations
		self.seed = 0

		#Which lensing quantities do we need?
		self.tomographic_convergence = False
		self.convergence = True
		self.shear = False
		self.omega = False

	def _init_plane_set(self):

		#Set of lens planes to be used during ray tracing
		self.plane_set = "Planes"
		self.plane_info_file = None

	def _init_randomizer(self):

		#N-body simulation realizations that need to be mixed
		self.mix_nbody_realizations = [1]
		self.mix_cut_points = [0]
		self.mix_normals = [0]
		self.lens_map_realizations = 4
		self.first_realization = 1

	###############################################################################################################################################

	@classmethod
	def get(cls,options):

		#Check that the config file has the appropriate section
		section = cls._section
		assert options.has_section(section),"No {0} section in configuration file {1}".format(section,options.filename)

		#Fill in the appropriate fields
		settings = cls()

		settings._read_common(options,section)
		settings._read_plane_set(options,section)
		settings._read_randomizer(options,section)

		#Return to user
		return settings

	def _read_common(self,options,section):

		self.directory_name = options.get(section,"directory_name")
		self.override_with_local = options.getboolean(section,"override_with_local")
		self.format = options.get(section,"format")
		
		try:
			self.plane_format = options.get(section,"plane_format")
		except NoOptionError:
			pass

		try:
			self.plane_name_format = options.get(section,"plane_name_format")
		except NoOptionError:
			pass
		
		self.map_resolution = options.getint(section,"map_resolution")
		
		self.angle_unit = getattr(u,options.get(section,"angle_unit"))
		self.map_angle = options.getfloat(section,"map_angle") * self.angle_unit
		
		self.source_redshift = options.getfloat(section,"source_redshift")

		self.seed = options.getint(section,"seed")

		try:
			self.tomographic_convergence = options.getboolean(section,"tomographic_convergence")
		except NoOptionError:
			pass

		self.convergence = options.getboolean(section,"convergence")
		self.shear = options.getboolean(section,"shear")
		self.omega = options.getboolean(section,"omega")

	def _read_plane_set(self,options,section):
		
		self.plane_set = options.get(section,"plane_set")
		try:
			self.plane_info_file = options.get(section,"plane_info_file")
		except NoOptionError:
			self.plane_info_file = None

	def _read_randomizer(self,options,section):

		self.mix_nbody_realizations = [ int(n) for n in options.get(section,"mix_nbody_realizations").split(",") ]
		self.lens_map_realizations = options.getint(section,"lens_map_realizations")
		self.mix_cut_points = [ int(n) for n in options.get(section,"mix_cut_points").split(",") ]
		self.mix_normals = [ int(n) for n in options.get(section,"mix_normals").split(",") ] 

		try:
			self.first_realization = options.getint(section,"first_realization")
		except NoOptionError:
			self.first_realization = 1



###########################################################
###########TelescopicMapSettings class#####################
###########################################################

class TelescopicMapSettings(MapSettings):

	"""
	Class handler of telescopic simulation map generation settings
	
	"""

	_section = "TelescopicMapSettings"

	def _init_plane_set(self):

		#Set of lens planes to be used during ray tracing
		self.plane_set = ("Planes",)
		self.plane_info_file = None

	def _init_randomizer(self):

		#N-body simulation realizations that need to be mixed
		self.mix_nbody_realizations = ([1],)
		self.mix_cut_points = ([0],)
		self.mix_normals = ([0],)
		self.lens_map_realizations = 4

	def _read_plane_set(self,options,section):
		
		self.plane_set = tuple(options.get(section,"plane_set").split(","))
		try:
			self.plane_info_file = options.get(section,"plane_info_file")
		except NoOptionError:
			self.plane_info_file = None

	def _read_randomizer(self,options,section):

		self.mix_nbody_realizations = ast.literal_eval(options.get(section,"mix_nbody_realizations"))
		self.lens_map_realizations = options.getint(section,"lens_map_realizations")
		self.mix_cut_points = ast.literal_eval(options.get(section,"mix_cut_points"))
		self.mix_normals = ast.literal_eval(options.get(section,"mix_normals"))

		try:
			self.first_realization = options.getint(section,"first_realization")
		except NoOptionError:
			self.first_realization = 1



#####################################################
###########CatalogSettings class#####################
#####################################################

class CatalogSettings(LTSettings):

	"""
	Class handler of simulated catalog generation settings

	"""

	def __init__(self,**kwargs):
	
		#Name of catalog batch
		self.directory_name = "Catalog"
		self.input_files = "galaxy_positions.fits"
		self.total_num_galaxies = 1000
		self.catalog_angle_unit = u.deg

		#Use the options generated at the moment of the batch generation (advised)
		self.override_with_local = True

		#Format of the simulated catalog files
		self.format = "fits"
		self.plane_format = "fits"
		self.plane_name_format = "snap{0}_potentialPlane{1}_normal{2}.{3}"

		#Random seed used to generate multiple catalog realizations
		self.seed = 0

		#Set of lens planes to be used during ray tracing
		self.plane_set = "Planes"

		#N-body simulation realizations that need to be mixed
		self.mix_nbody_realizations = [1]
		self.mix_cut_points = [0]
		self.mix_normals = [0]
		self.lens_catalog_realizations = 1
		self.first_realization = 1

		#Allow for kwargs override
		for key in kwargs.keys():
			setattr(self,key,kwargs[key])


	@classmethod
	def get(cls,options):

		#Check that the config file has the appropriate section
		section = "CatalogSettings"
		assert options.has_section(section),"No {0} section in configuration file {1}".format(section,options.filename)

		#Fill in the appropriate fields
		settings = cls()

		#Name of catalog batch
		settings.directory_name = options.get(section,"directory_name")
		settings.input_files = options.get(section,"input_files").split(",")
		settings.total_num_galaxies = options.getint(section,"total_num_galaxies")
		settings.catalog_angle_unit = getattr(u,options.get(section,"catalog_angle_unit"))

		#Use the options generated at the moment of the batch generation (advised)
		settings.override_with_local = options.getboolean(section,"override_with_local")

		#Format of the simulated catalog files
		settings.format = options.get(section,"format")

		try:
			self.plane_format = options.get(section,"plane_format")
		except NoOptionError:
			pass

		try:
			self.plane_name_format = options.get(section,"plane_name_format")
		except NoOptionError:
			pass

		#Set of lens planes to be used during ray tracing
		settings.seed = options.getint(section,"seed")

		#Set of lens planes to be used during ray tracing
		settings.plane_set = options.get(section,"plane_set")

		#N-body simulation realizations that need to be mixed
		settings.mix_nbody_realizations = [ int(n) for n in options.get(section,"mix_nbody_realizations").split(",") ]
		settings.lens_catalog_realizations = options.getint(section,"lens_catalog_realizations")
		
		try:
			settings.realizations_per_subdirectory = options.getint(section,"realizations_per_subdirectory")
		except NoOptionError:
			pass
		
		settings.mix_cut_points = [ int(n) for n in options.get(section,"mix_cut_points").split(",") ]
		settings.mix_normals = [ int(n) for n in options.get(section,"mix_normals").split(",") ]

		try:
			settings.first_realization = options.getint(section,"first_realization")
		except NoOptionError:
			settings.first_realization = 1

		#Return to user
		return settings


##################################################
###############JobSettings class##################
##################################################

class JobSettings(LTSettings):

	"""
	Class handler of batch job submission settings

	"""

	def __init__(self,**kwargs):

		#Personal settings
		self.email = "apetri@phys.columbia.edu"
		self.charge_account = "TG-AST140041"

		#Path to executable
		self.path_to_executable = "Gadget2"

		#Name of the job, output
		self.job_name = "job"
		self.redirect_stdout = "job.out"
		self.redirect_stderr = "job.err"

		#Resources
		self.cores_per_simulation = 16
		self.queue = "development"
		self.wallclock_time = "02:00:00"

		#Script name
		self.job_script_file = "job.sh"

		#Allow for kwargs override
		for key in kwargs.keys():
			setattr(self,key,kwargs[key])


	@classmethod
	def get(cls,options,section):

		#Check that the config file has the appropriate section
		assert options.has_section(section),"No {0} section in configuration file {1}".format(section,options.filename)

		#Fill in the appropriate fields
		settings = cls()

		#Personal settings
		settings.email = options.get(section,"email")
		settings.charge_account = options.get(section,"charge_account")

		#Path to executable
		try:
			settings.path_to_executable = options.get(section,"path_to_executable")
		except NoOptionError:
			settings.path_to_executable = section

		#Name of the job, output
		settings.job_name = options.get(section,"job_name")
		settings.redirect_stdout = options.get(section,"redirect_stdout")
		settings.redirect_stderr = options.get(section,"redirect_stderr")

		#Resources
		
		#These do not need to be provided necessarily
		try:
			settings.num_cores = options.getint(section,"num_cores")
		except NoOptionError:
			pass

		try:
			settings.num_nodes = options.getint(section,"num_nodes")
		except NoOptionError:
			pass

		try:
			settings.tasks_per_node = options.getint(section,"tasks_per_node")
		except NoOptionError:
			pass

		#These need to be provided
		settings.cores_per_simulation = options.getint(section,"cores_per_simulation")
		settings.queue = options.get(section,"queue")
		settings.wallclock_time = options.get(section,"wallclock_time")

		#Script name
		settings.job_script_file = options.get(section,"job_script_file")

		return settings















