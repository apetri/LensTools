from __future__ import division
import sys,os

if sys.version_info.major>=3:
	from io import StringIO
else:
	from StringIO import StringIO

from .nbody import NbodySnapshot
from .. import extern as ext
from .settings import LTSettings

import numpy as np
import astropy.units as u
from astropy.cosmology import w0waCDM,LambdaCDM

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

		output = StringIO()

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

############################################################
################Gadget2Header class#########################
############################################################

class Gadget2Header(dict):

	"""
	Class handler of a Gadget2 snapshot header

	"""

	def __init__(self,HeaderDict=dict()):

		super(Gadget2Header,self).__init__()
		for key in HeaderDict.keys():
			self[key] = HeaderDict[key]

	def __repr__(self):

		keys = self.keys()
		keys.sort()
		
		return "\n".join([ "{0} : {1}".format(key,self[key]) for key in keys ]) 

	def __add__(self,rhs):

		assert isinstance(rhs,Gadget2Header),"addition not defined if rhs is not a Gadget2Header!"

		#Check that it makes sense to add the snapshots (cosmological parameters, box size, time and redshift must agree)
		fields_to_match = ["Ode0","Om0","h","w0","wa","box_size","endianness","flag_cooling","flag_feedback","flag_sfr","num_files"]
		fields_to_match += ["num_particles_total","num_particles_total_gas","num_particles_total_side","num_particles_total_with_mass","redshift","scale_factor"]

		for field in fields_to_match:
			assert self[field] == rhs[field],"{0} fields do not match!".format(field)

		assert np.all(self["masses"]==rhs["masses"])
		assert np.all(self["num_particles_total_of_type"]==rhs["num_particles_total_of_type"])

		#Construct the header of the merged snapshot
		merged_header = self.copy()
		merged_header["files"] += rhs["files"]
		merged_header["num_particles_file"] += rhs["num_particles_file"]
		merged_header["num_particles_file_gas"] += rhs["num_particles_file_gas"]
		merged_header["num_particles_file_of_type"] += rhs["num_particles_file_of_type"]
		merged_header["num_particles_file_with_mass"] += rhs["num_particles_file_with_mass"]

		return merged_header

##############################################################
#################Gadget2Snapshot class######################
##############################################################

class Gadget2Snapshot(NbodySnapshot):

	"""
	A class that handles Gadget2 snapshots, mainly I/O from the binary format and spatial information statistics.Inherits from the abstract NbodySnapshot

	"""

	###############################################################################################
	#########################Abstract methods implementation#######################################
	###############################################################################################

	@classmethod
	def buildFilename(cls,root,pool,**kwargs):
		
		if pool is not None:
			return root+".{0}".format(pool.rank)
		else:
			return root

	@classmethod
	def int2root(cls,name,n):
		return name + "_{0:03d}".format(n)

	############################################################################################

	def getHeader(self):
		
		header = Gadget2Header(ext._gadget2.getHeader(self.fp))
		header["files"] = [self.fp.name]

		header["w0"] = -1.0
		header["wa"] = 0.0
		header["comoving_distance"] = LambdaCDM(H0=header["h"]*100,Om0=header["Om0"],Ode0=header["Ode0"]).comoving_distance(header["redshift"]).to(u.kpc).value * header["h"]

		return header 

	############################################################################################

	def setLimits(self):
		self._first = None
		self._last = None

	############################################################################################

	def getPositions(self,first=None,last=None,save=True):

		"""
		Reads in the particles positions (read in of a subset is allowed): when first and last are specified, the numpy array convention is followed (i.e. getPositions(first=a,last=b)=getPositions()[a:b])

		:param first: first particle in the file to be read, if None 0 is assumed
		:type first: int. or None

		:param last: last particle in the file to be read, if None the total number of particles is assumed
		:type last: int. or None

		:param save: if True saves the particles positions as attribute
		:type save: bool.

		:returns: numpy array with the particle positions

		"""

		#The file must not be closed
		assert not self.fp.closed

		#Particles do not have structure
		self.weights = None
		self.virial_radius = None
		self.concentration = None

		numPart = self._header["num_particles_file"]

		#Calculate the offset from the beginning of the file: 4 bytes (endianness) + 256 bytes (header) + 8 bytes (void)
		offset = 4 + 256 + 8

		#If first is specified, offset the file pointer by that amount
		if first is not None:
			
			assert first>=0
			offset += 4 * 3 * first
			numPart -= first

		if last is not None:

			if first is not None:
				
				assert last>=first and last<=self._header["num_particles_file"]
				numPart = last - first

			else:

				assert last<=self._header["num_particles_file"]
				numPart = last


		#Read in the particles positions and return the corresponding array
		try:
			positions = (ext._gadget2.getPosVel(self.fp,offset,numPart) * self.kpc_over_h).to(self.Mpc_over_h)
		except AttributeError:
			positions = ext._gadget2.getPosVel(self.fp,offset,numPart) * u.kpc

		if save:
			self.positions = positions
			return self.positions
		
		#Return
		return positions

	############################################################################################

	def getVelocities(self,first=None,last=None,save=True):

		"""
		Reads in the particles velocities (read in of a subset is allowed): when first and last are specified, the numpy array convention is followed (i.e. getVelocities(first=a,last=b)=getVelocities()[a:b])

		:param first: first particle in the file to be read, if None 0 is assumed
		:type first: int. or None

		:param last: last particle in the file to be read, if None the total number of particles is assumed
		:type last: int. or None

		:param save: if True saves the particles velocities as attrubute
		:type save: bool.

		:returns: numpy array with the particle velocities

		"""

		assert not self.fp.closed

		numPart = self._header["num_particles_file"]

		#Calculate the offset from the beginning of the file: 4 bytes (endianness) + 256 bytes (header) + 8 bytes (void)
		offset = 4 + 256 + 8

		#Skip all the particle positions
		offset += 4 * 3 * numPart

		#Skip other 8 void bytes
		offset += 8

		#If first is specified, offset the file pointer by that amount
		if first is not None:
			
			assert first>=0
			offset += 4 * 3 * first
			numPart -= first

		if last is not None:

			if first is not None:
				
				assert last>=first and last<=self._header["num_particles_file"]
				numPart = last - first

			else:

				assert last<=self._header["num_particles_file"]
				numPart = last


		#Read in the particles positions and return the corresponding array
		velocities = ext._gadget2.getPosVel(self.fp,offset,numPart)

		#Scale units
		velocities *= self._velocity_unit
		velocities *= u.cm / u.s

		if save:
			self.velocities = velocities
			return self.velocities
		
		#Return
		return velocities

	############################################################################################

	def getID(self,first=None,last=None,save=True):

		"""
		Reads in the particles IDs, 4 byte ints, (read in of a subset is allowed): when first and last are specified, the numpy array convention is followed (i.e. getID(first=a,last=b)=getID()[a:b])

		:param first: first particle in the file to be read, if None 0 is assumed
		:type first: int. or None

		:param last: last particle in the file to be read, if None the total number of particles is assumed
		:type last: int. or None

		:param save: if True saves the particles IDs as attribute
		:type save: bool.

		:returns: numpy array with the particle IDs

		"""

		assert not self.fp.closed

		numPart = self._header["num_particles_file"]

		#Calculate the offset from the beginning of the file: 4 bytes (endianness) + 256 bytes (header) + 8 bytes (void)
		offset = 4 + 256 + 8

		#Skip all the particle positions
		offset += 4 * 3 * numPart

		#Skip other 8 void bytes
		offset += 8

		#Skip all the particle velocities
		offset += 4 * 3 * numPart

		#Skip other 8 void bytes
		offset += 8

		#If first is specified, offset the file pointer by that amount
		if first is not None:
			
			assert first>=0
			offset += 4 * first
			numPart -= first

		if last is not None:

			if first is not None:
				
				assert last>=first and last<=self._header["num_particles_file"]
				numPart = last - first

			else:

				assert last<=self._header["num_particles_file"]
				numPart = last


		#Read in the particles positions and return the corresponding array
		ids = ext._gadget2.getID(self.fp,offset,numPart)
		if save:
			self.id = ids
			return self.id
		
		#Return
		return ids

	############################################################################################

	def write(self,filename,files=1):

		"""
		Writes particles information (positions, velocities, etc...) to a properly formatter Gadget snapshot

		:param filename: name of the file to which to write the snapshot
		:type filename: str.

		:param files: number of files on which to split the writing of the snapshot (useful if the number of particles is large); if > 1 the extension ".n" is appended to the filename
		:type files: int.

		"""

		#Sanity checks
		assert hasattr(self,"positions"),"Positions must be specified!!"
		assert self.positions.shape[1]==3

		if not hasattr(self,"_header"):
			self.setHeaderInfo()	

		#Build a bare header based on the available info (need to convert units back to the Gadget ones)
		_header_bare = self._header.copy()
		_header_bare["box_size"] = _header_bare["box_size"].to(self.kpc_over_h).value
		_header_bare["masses"] = _header_bare["masses"].to(u.g).value * _header_bare["h"] / self._mass_unit
		_header_bare["num_particles_file_of_type"] = _header_bare["num_particles_file_of_type"].astype(np.int32)
		_header_bare["num_particles_total_of_type"] = _header_bare["num_particles_total_of_type"].astype(np.int32)
		_header_bare["comoving_distance"] = _header_bare["comoving_distance"].to(self.Mpc_over_h).value * 1.0e3


		#Convert units for positions and velocities
		_positions_converted = self.positions.to(self.kpc_over_h).value.astype(np.float32)

		if hasattr(self,"velocities"):
			assert self.positions.shape==self.velocities.shape
			_velocities_converted = (self.velocities.to(u.cm/u.s).value / self._velocity_unit).astype(np.float32)
			writeVel = 1
		else:
			_velocities_converted = np.zeros((1,3),dtype=np.float32)
			writeVel = 0

		#Check if we want to split on multiple files (only DM particles supported so far for this feature)
		if files>1:

			#Number of files the snapshot is split into
			_header_bare["num_files"] = files

			#Update the header with the file names
			self.header["files"] = [ "{0}.{1}".format(filename,n) for n in range(files) ]

			#Distribute particles among files
			particles_per_file = _header_bare["num_particles_total"] // files

			#Write each file
			for n in range(files - 1):
				
				#Update header
				_header_bare["num_particles_file"] = particles_per_file
				#TODO all particles are DM, fix distribution in the future
				_header_bare["num_particles_file_of_type"] = np.array([0,particles_per_file,0,0,0,0],dtype=np.int32)

				#Write it!
				filename_with_extension = "{0}.{1}".format(filename,n)
				ext._gadget2.write(_header_bare,_positions_converted[n*particles_per_file:(n+1)*particles_per_file],_velocities_converted[n*particles_per_file:(n+1)*particles_per_file],n*particles_per_file+1,filename_with_extension,writeVel)

			
			#The last file might have a different number of particles
			particles_last_file = len(_positions_converted[particles_per_file*(files-1):])

			#Update header
			_header_bare["num_particles_file"] = particles_last_file
			#TODO all particles are DM, fix distribution in the future
			_header_bare["num_particles_file_of_type"] = np.array([0,particles_last_file,0,0,0,0],dtype=np.int32)

			#Write it!
			filename_with_extension = "{0}.{1}".format(filename,files-1)
			ext._gadget2.write(_header_bare,_positions_converted[particles_per_file*(files-1):],_velocities_converted[particles_per_file*(files-1):],(files-1)*particles_per_file+1,filename_with_extension,writeVel)

		else:

			#Update the num_files key only if not present already
			if "num_files" not in _header_bare.keys():
				_header_bare["num_files"] = 1

			#Update the header with the file names
			self.header["files"] = [ filename ]
			
			#Write it!!
			ext._gadget2.write(_header_bare,_positions_converted,_velocities_converted,1,filename,writeVel)

	############################################################################################
	###########################Extra methods####################################################
	############################################################################################

	def setHeaderInfo(self,Om0=0.26,Ode0=0.74,w0=-1.0,wa=0.0,h=0.72,redshift=100.0,box_size=15.0*u.Mpc/0.72,flag_cooling=0,flag_sfr=0,flag_feedback=0,flag_stellarage=0,flag_metals=0,flag_entropy_instead_u=0,masses=np.array([0,1.03e10,0,0,0,0])*u.Msun,num_particles_file_of_type=None,npartTotalHighWord=np.zeros(6,dtype=np.uint32)):

		"""
		Sets the header info in the snapshot to write

		"""

		if num_particles_file_of_type is None:
			num_particles_file_of_type = np.array([0,1,0,0,0,0],dtype=np.int32) * self.positions.shape[0]

		assert num_particles_file_of_type.sum()==self.positions.shape[0],"The total number of particles must match!!"
		assert box_size.unit.physical_type=="length"
		assert masses.unit.physical_type=="mass"

		#Create the header
		self._header = Gadget2Header()
		
		#Fill in
		self._header["Om0"] = Om0
		self._header["Ode0"] = Ode0
		self._header["w0"] = w0
		self._header["wa"] = wa
		self._header["h"] = h
		self._header["H0"] = 100.0*h*u.km/(u.s*u.Mpc)
		self._header["redshift"] = redshift
		self._header["scale_factor"] = 1.0 / (1.0 + redshift)
		self._header["box_size"] = box_size
		self._header["flag_cooling"] = flag_cooling
		self._header["flag_sfr"] = flag_sfr
		self._header["flag_feedback"] = flag_feedback
		self._header["flag_stellarage"] = flag_stellarage
		self._header["flag_metals"] = flag_metals
		self._header["flag_entropy_instead_u"] = flag_entropy_instead_u
		self._header["masses"] = masses
		self._header["num_particles_file_of_type"] = num_particles_file_of_type
		self._header["num_particles_file"] = num_particles_file_of_type.sum()
		self._header["num_particles_total_of_type"] = num_particles_file_of_type
		self._header["num_particles_total"] = num_particles_file_of_type.sum()
		self._header["npartTotalHighWord"] = npartTotalHighWord

		#Define the kpc/h and Mpc/h units for convenience
		self.kpc_over_h = u.def_unit("kpc/h",u.kpc/self._header["h"])
		self.Mpc_over_h = u.def_unit("Mpc/h",u.Mpc/self._header["h"])

		#Compute the comoving distance according to the model
		cosmo = w0waCDM(H0=100.0*h,Om0=Om0,Ode0=Ode0,w0=w0,wa=wa)
		self._header["comoving_distance"] = cosmo.comoving_distance(redshift).to(self.Mpc_over_h)


	def writeParameterFile(self,filename,settings):

		"""
		Writes a Gadget2 parameter file to evolve the current snapshot using Gadget2

		:param filename: name of the file to which to write the parameters
		:type filename: str.

		:param settings: tunable settings of Gadget2 (see Gadget2 manual)
		:type settings: Gadget2Settings

		"""

		#Create output directory if not existent already
		outputdir = settings.OutputDir
		if not(os.path.isdir(outputdir)):
			os.mkdir(outputdir)

		#Update the settings according to the physical units of the current snapshot
		settings.UnitLength_in_cm = self._length_unit
		settings.UnitMass_in_g = self._mass_unit
		settings.UnitVelocity_in_cm_per_s = self._velocity_unit

		#Set the appropriate name for the initial condition file
		if "files" in self.header.keys():
			suffix = self.header["files"][0].split(".")[-1]
			settings.InitCondFile = os.path.abspath(self.header["files"][0].rstrip(".{0}".format(suffix)))

		#Write the options
		with open(filename,"w") as paramfile:

			#Filenames section
			paramfile.write(settings.writeSection("file_names"))
			
			#CPU time limit section
			paramfile.write(settings.writeSection("cpu_timings"))

			#Code options section
			paramfile.write(settings.writeSection("code_options"))

			#Initial scale factor time
			paramfile.write("{0}			{1}\n".format("TimeBegin",self.header["scale_factor"]))

			#Characteristics of run section
			paramfile.write(settings.writeSection("characteristics_of_run"))
			
			#Cosmological parameters
			paramfile.write("{0}			{1}\n".format("Omega0",self.header["Om0"]))
			paramfile.write("{0}			{1}\n".format("OmegaLambda",self.header["Ode0"]))
			paramfile.write("{0}			{1}\n".format("OmegaBaryon",0.046))
			paramfile.write("{0}			{1}\n".format("HubbleParam",self.header["h"]))
			paramfile.write("{0}			{1}\n".format("BoxSize",self.header["box_size"].to(self.kpc_over_h).value))
			paramfile.write("{0}			{1}\n".format("w0",self.header["w0"]))
			paramfile.write("{0}			{1}\n\n".format("wa",self.header["wa"]))

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


##############################################################
#################Gadget2SnapshotDE class########################
##############################################################

class Gadget2SnapshotDE(Gadget2Snapshot):

	"""
	A class that handles Gadget2 snapshots, mainly I/O from the binary format and spatial information statistics.Inherits from Gadget2Snapshot; assumes that the header includes Dark Energy information

	"""

	def getHeader(self):
		header = Gadget2Header(ext._gadget2.getHeader(self.fp))
		header["files"] = [self.fp.name]

		return header
		 

##################################################################
#################Gadget2SnapshotPipe class########################
##################################################################

class Gadget2SnapshotPipe(Gadget2SnapshotDE):

	"""
	Read in the particle positions when calling the constructor, without calling fseek

	"""

	def __init__(self,*args,**kwargs):

		#Call parent constructor
		super(Gadget2SnapshotPipe,self).__init__(*args,**kwargs)

		#Read in the positions 
		npart = self.header["num_particles_file"]
		self.fp.read(8)
		
		try:
			self.positions = (np.fromstring(self.fp.read(4*3*npart),dtype=np.float32).reshape(npart,3) * self.kpc_over_h).to(self.Mpc_over_h)
		except AttributeError:
			pass

		#Read the rest
		self.fp.read()

		#Particles do not have structure
		self.weights = None
		self.virial_radius = None
		self.concentration = None


