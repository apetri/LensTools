import sys,os
import ast

if sys.version_info.major>=3:
	from configparser import NoOptionError
else:
	from ConfigParser import NoOptionError

import json

import numpy as np
import astropy.units as u

from ..simulations.settings import select_parser,LTSettings
from ..simulations.camb import CAMBSettings
from ..simulations.gadget2 import Gadget2Settings


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
		for key in kwargs:
			setattr(self,key,kwargs[key])


###########################################################################################################################
###########################################################################################################################			

#################################################
###########PlaneSettings class###################
#################################################

class PlaneSettings(LTSettings):

	"""
	Class handler of plane generation settings from constant time Nbody snapshots

	"""

	_section = "PlaneSettings"

	def __init__(self,**kwargs):

		#Name of the planes batch
		self.directory_name = "Planes"

		#Snapshot class handler
		self.snapshot_handler = "Gadget2SnapshotPipe"
		
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
		for key in kwargs:
			setattr(self,key,kwargs[key])

	@classmethod
	def get(cls,options):

		#Check that the config file has the appropriate section
		section = cls._section
		assert options.has_section(section),"No {0} section in configuration file {1}".format(section,options.filename)

		#Fill in the appropriate fields
		settings = cls()
		
		settings.directory_name = options.get(section,"directory_name")

		#Snapshot class handler
		try:
			settings.snapshot_handler = options.get(section,"snapshot_handler")
		except NoOptionError:
			pass

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

##########################################################
###########PlaneLightConeSettings class###################
##########################################################

class PlaneLightConeSettings(LTSettings):

	"""
	Class handler of plane generation settings from lightcone projection Nbody snapshots

	"""

	_section = "PlaneLightConeSettings"

	def __init__(self,**kwargs):

		#Name of the planes batch
		self.directory_name = "Planes"

		#Snapshot class handler
		self.snapshot_handler = "FastPMSnapshot"
		
		#Use the pickled options generated at the moment of the batch generation (advised)
		self.override_with_local = False

		self.format = "fits"
		self.name_format = "snap{0}_{1}Plane{2}_normal{3}.{4}"

		#Lens discretization
		self.zmax = 3.0
		self.num_lenses = 20
		self.normal = 2
		self.plane_resolution = 64

		#Optional, not usually changed
		self.thickness_resolution = 1
		self.smooth = 1
		self.kind = "potential"

		#On the fly raytracing
		self.do_lensing = False
		self.integration_type = "full"
		self.fov = 3.0*u.deg
		self.fov_resolution = 32

		#Allow for kwargs override
		for key in kwargs:
			setattr(self,key,kwargs[key])

	@classmethod
	def get(cls,options):

		#Check that the config file has the appropriate section
		section = cls._section
		assert options.has_section(section),"No {0} section in configuration file {1}".format(section,options.filename)

		#Fill in the appropriate fields
		settings = cls()
		
		settings.directory_name = options.get(section,"directory_name")

		#Snapshot class handler
		try:
			settings.snapshot_handler = options.get(section,"snapshot_handler")
		except NoOptionError:
			pass

		settings.override_with_local = options.getboolean(section,"override_with_local")
		
		settings.format = options.get(section,"format")
		try:
			settings.name_format = options.get(section,"name_format")
		except NoOptionError:
			pass
		
		settings.plane_resolution = options.getint(section,"plane_resolution")

		#Lens discretization
		settings.zmax = options.getfloat(section,"zmax")
		settings.num_lenses = options.getint(section,"num_lenses")
		settings.normal = options.getint(section,"normal")
		settings.plane_resolution = options.getint(section,"plane_resolution")

		#Weak lensing
		settings.do_lensing = options.getboolean(section,"do_lensing")
		settings.integration_type = options.get(section,"integration_type")
		settings.fov = options.getfloat(section,"fov_deg")*u.deg
		settings.fov_resolution = options.getint(section,"fov_resolution")

		###########
		#Optionals#
		###########

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


###########################################################################################################################
###########################################################################################################################	

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
		for key in kwargs:
			setattr(self,key,kwargs[key])

	def _init_commmon(self):

		#Names of the map batch
		self.directory_name = "Maps"

		#Use the options generated at the moment of the batch generation (advised)
		self.override_with_local = True

		self.format = "fits"
		self.plane_format = "fits"
		self.plane_name_format = "snap{0}_potentialPlane{1}_normal{2}.{3}"
		self.lens_type = "PotentialPlane"

		self.map_resolution = 128
		self.map_angle = 1.6 * u.deg
		self.angle_unit = u.deg
		self.source_redshift = 2.0

		#Random seed used to generate multiple map realizations
		self.seed = 0

		#Transpose lenses up to a certain index
		self.transpose_up_to = -1

		#Which lensing quantities do we need?
		self.tomographic_convergence = False
		self.convergence = True
		self.convergence_ks = False
		self.shear = False
		self.reduced_shear = False
		self.reduced_shear_convergence = False
		self.omega = False

		#Line of sight integration type
		self.integration_type = "born"

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

		try:
			self.lens_type = options.get(section,"lens_type")
		except NoOptionError:
			pass
		
		self.map_resolution = options.getint(section,"map_resolution")
		
		self.angle_unit = getattr(u,options.get(section,"angle_unit"))
		self.map_angle = options.getfloat(section,"map_angle") * self.angle_unit
		
		self.source_redshift = options.getfloat(section,"source_redshift")

		self.seed = options.getint(section,"seed")
		try:
			self.transpose_up_to = options.getint(section,"transpose_up_to")
		except NoOptionError:
			pass

		###########################################################################################

		try:
			self.tomographic_convergence = options.getboolean(section,"tomographic_convergence")
		except NoOptionError:
			pass

		try:
			self.convergence = options.getboolean(section,"convergence")
		except NoOptionError:
			pass

		try:
			self.convergence_ks = options.getboolean(section,"convergence_ks")
		except NoOptionError:
			pass

		try:
			self.shear = options.getboolean(section,"shear")
		except NoOptionError:
			pass

		try:
			self.reduced_shear = options.getboolean(section,"reduced_shear")
		except NoOptionError:
			pass

		try:
			self.reduced_shear_convergence = options.getboolean(section,"reduced_shear_convergence")
		except NoOptionError:
			pass

		try:
			self.omega = options.getboolean(section,"omega")
		except NoOptionError:
			pass

		###########################################################################################

		try:
			self.integration_type = options.get(section,"integration_type")
		except NoOptionError:
			pass

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

	_section = "CatalogSettings"

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

		#Reduced shear
		self.reduced_shear = False

		#Allow for kwargs override
		for key in kwargs:
			setattr(self,key,kwargs[key])


	@classmethod
	def get(cls,options):

		#Check that the config file has the appropriate section
		section = cls._section
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

		try:
			settings.reduced_shear = options.getboolean(section,"reduced_shear")
		except NoOptionError:
			pass

		#Return to user
		return settings

#####################################################
###########CMBReconstructionSettings class###########
#####################################################

class CMBReconstructionSettings(LTSettings):

	"""
	Class handler of CMB lensing reconstruction settings

	"""

	_section = "CMBReconstruction"

	def __init__(self,**kwargs):

		#Input
		self.input_set = "kappaCMB"
		self.input_filename = "WLconv*.fits"
		
		#Quadratic estimator settings
		self.estimator = "TT"
		self.unlensed_ps_filename = "scals.dat"
		self.lensed_ps_filename = "lensed_scals.dat"
		self.ps_type = "camb_dimensionless"
		self.lmax = 3500

		#What to estimate
		self.output_type = "kappa"

		#Noise/filtering
		self.wiener = False
		self.noise_level = 6.0*u.uK*u.arcmin
		self.beam_fwhm = 1.4*u.arcmin
		
		#Output
		self.save_intermediate = False
		self.output_set = "kappaCMBRec"
		self.output_fname = "{0}_{1}r.fits"

		#Allow for kwargs override
		for key in kwargs:
			setattr(self,key,kwargs[key])

	@classmethod
	def get(cls,options):

		#Check that the config file has the appropriate section
		section = cls._section
		assert options.has_section(section),"No {0} section in configuration file {1}".format(section,options.filename)

		#Fill in the appropriate fields
		settings = cls()

		#Input
		settings.input_set = options.get(section,"input_set")
		settings.input_filename = options.get(section,"input_filename")
		
		#Quadratic estimator settings
		settings.estimator = options.get(section,"estimator")
		settings.unlensed_ps_filename = options.get(section,"unlensed_ps_filename")
		settings.lensed_ps_filename = options.get(section,"lensed_ps_filename")
		settings.ps_type = options.get(section,"ps_type")
		settings.lmax = options.getint(section,"lmax")

		#What to estimate
		settings.output_type = options.get(section,"output_type")

		#Noise/filtering
		settings.wiener = options.getboolean(section,"wiener")
		settings.noise_level = options.getfloat(section,"noise_level_uK_arcmin")*u.uK*u.arcmin
		settings.beam_fwhm = options.getfloat(section,"beam_fwhm_arcmin")*u.arcmin
		
		#Output
		settings.save_intermediate = options.getboolean(section,"save_intermediate")
		settings.output_set = options.get(section,"output_set")
		settings.output_fname = options.get(section,"output_fname")

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
		for key in kwargs:
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















