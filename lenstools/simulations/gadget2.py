from __future__ import division
from operator import mul
from functools import reduce

import os
import StringIO

from .. import extern as ext

import numpy as np

#astropy stuff, invaluable here
from astropy.units import Mbyte,kpc,Mpc,cm,km,g,s,hour,day,deg,arcmin,rad,Msun,quantity,def_unit
from astropy.constants import c
from astropy.cosmology import w0waCDM

#FFT engines
from numpy.fft import rfftn,irfftn,fftfreq
try:
	from numpy.fft import rfftfreq
except ImportError:
	from .. import utils
	rfftfreq = utils.rfftfreq

#KD-Tree
from scipy.spatial import cKDTree as KDTree

#Plotting engine
try:
	import matplotlib.pyplot as plt
	from mpl_toolkits.mplot3d import Axes3D
	matplotlib = True
except ImportError:

	matplotlib = False

#Try to import r2py to save snapshot positions in R format
try:
	import rpy2.robjects as robj
	rpy2 = True
except ImportError:
	rpy2 = False

############################################################
################Gadget2Settings class#######################
############################################################

class Gadget2Settings(object):

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

	def __init__(self):

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
		self.TimeLimitCPU = 1.0*day
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
		self.CpuTimeBetRestartFile = 12.5*hour 
		self.TimeBetStatistics = 0.05
		self.NumFilesPerSnapshot = 1
		self.NumFilesWrittenInParallel = 1

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
		self.PartAllocFactor = 2    
		self.TreeAllocFactor = 1 
		self.BufferSize = 20*Mbyte 


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
			if type(value)==quantity.Quantity:
				
				if value.unit.physical_type=="time":
					value = value.to(s).value
				elif value.unit.physical_type=="speed":
					value = value.to(cm/s).value
				elif "byte" in value.unit.to_string():
					value = value.to(Mbyte).value

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

############################################################
#################Gadget2Snapshot class######################
############################################################

class Gadget2Snapshot(object):

	"""
	A class that handles Gadget2 snapshots, mainly I/O from the binary format and spatial information statistics

	"""

	def __init__(self,fp=None,pool=None,length_unit=1.0*kpc,mass_unit=1.0e10*Msun,velocity_unit=1.0*km/s):

		self._length_unit = length_unit.to(cm).value
		self._mass_unit = mass_unit.to(g).value
		self._velocity_unit = velocity_unit.to(cm/s).value

		assert (type(fp)==file) or (fp is None),"Call the open() method instead!!"

		if fp is not None:
		
			self.fp = fp

			self._header = Gadget2Header(ext._gadget.getHeader(fp))
			self._header["files"] = [self.fp.name]
			h = self._header["h"]

			#Define the Mpc/h, and kpc/h units for convenience
			if h>0.0:
				
				self.kpc_over_h = def_unit("kpc/h",kpc/self._header["h"])
				self.Mpc_over_h = def_unit("Mpc/h",Mpc/self._header["h"])

				#Scale box to kpc/h
				self._header["box_size"] *= self.kpc_over_h
				#Convert to Mpc/h
				self._header["box_size"] = self._header["box_size"].to(self.Mpc_over_h)

				#Read in the comoving distance
				self._header["comoving_distance"] = (self._header["comoving_distance"] / 1.0e3) * self.Mpc_over_h

			else:
				self._header["box_size"] *= kpc
				print("Warning! Hubble parameter h is zero!!")

			#Scale masses to correct units
			if h>0.0:
				self._header["masses"] *= (self._mass_unit / self._header["h"])
				self._header["masses"] *= g
				self._header["masses"] = self._header["masses"].to(Msun) 

			#Scale Hubble parameter to correct units
			self._header["H0"] = self._header["h"] * 100 * km / (s*Mpc)

			#Update the dictionary with the number of particles per side
			self._header["num_particles_total_side"] = int(np.round(self._header["num_particles_total"]**(1/3)))

			#Once all the info is available, add a wCDM instance as attribute to facilitate the cosmological calculations
			if h>0.0:
				self.cosmology = w0waCDM(H0=self._header["H0"],Om0=self._header["Om0"],Ode0=self._header["Ode0"],w0=self._header["w0"],wa=self._header["wa"])

		self.pool = pool

	@classmethod
	def open(cls,filename,pool=None):

		"""
		Opens a gadget snapshot at filename

		:param filename: file name of the gadget snapshot
		:type filename: str. or file.

		:param pool: use to distribute the calculations on different processors; if not None, each processor takes care of one of the snapshot parts, appending as ".n" to the filename
		:type pool: MPIWhirlPool instance

		"""

		if type(filename)==str:

			if pool is not None:
				filename+=".{0}".format(pool.rank)
			
			fp = open(filename,"r")
		
		elif type(filename)==file:
			
			if pool is not None:
				raise TypeError("Specifying file objects with MPIPools is not allowed!")
			fp = filename
		
		else:
			raise TypeError("filename must be string or file!")
		
		return cls(fp,pool)

	@property
	def header(self):

		"""
		Displays the snapshot header information

		:returns: the snapshot header information in dictionary form
		:rtype: dict.

		"""

		return self._header

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

		assert not self.fp.closed

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
			positions = (ext._gadget.getPosVel(self.fp,offset,numPart) * self.kpc_over_h).to(self.Mpc_over_h)
		except AttributeError:
			positions = ext._gadget.getPosVel(self.fp,offset,numPart) * kpc

		if save:
			self.positions = positions
		
		#Return
		return positions


	def pos2R(self,filename,variable_name="pos"):

		"""
		Saves the positions of the particles in a R readable format, for facilitating visualization with RGL

		:param filename: name of the file on which to save the particles positions
		:type filename: str.

		:param variable_name: name of the variable that contains the (x,y,z) positions in the R environment
		:type variable_name: str.

		"""

		if not rpy2:
			raise ImportError("rpy2 is not installed, can't proceed!")

		#Read in the positions
		if not hasattr(self,"positions"):
			self.getPositions()

		#Convert numpy array into an R vector
		positions_bare = self.positions.to(Mpc).value
		r_positions = robj.FloatVector(positions_bare.T.ravel())

		#Set the R environment
		robj.rinterface.globalenv[variable_name] = robj.r["matrix"](r_positions,nrow=positions_bare.shape[0])

		#Save
		robj.r.save(variable_name,file=filename)

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
		velocities = ext._gadget.getPosVel(self.fp,offset,numPart)

		#Scale units
		velocities *= self._velocity_unit
		velocities *= cm / s

		if save:
			self.velocities = velocities
		
		#Return
		return velocities

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
		ids = ext._gadget.getID(self.fp,offset,numPart)
		if save:
			self.id = ids
		
		#Return
		return ids


	def reorder(self):

		"""
		Sort particles attributes according to their ID

		"""

		assert hasattr(self,"id")
		
		#Rank the IDs
		idx = np.argsort(self.id)

		#Sort positions
		if hasattr(self,"positions"):
			
			assert self.positions.shape[0]==len(self.id)
			self.positions = self.positions[idx]

		#Sort velocities
		if hasattr(self,"velocities"):

			assert self.velocities.shape[0]==len(self.id)
			self.velocities = self.velocities[idx]

		#Finally sort IDs
		self.id.sort()


	def gridID(self):

		"""
		Compute an ID for the particles in incresing order according to their position on a Nside x Nside x Nside grid; the id is computed as x + y*Nside + z*Nside**2

		:returns: the gridded IDs
		:rtype: array of float

		"""

		try:
			pos = self.positions
		except:
			pos = self.getPositions()

		#Set the measure units for the grid
		grid_unit = self.header["box_size"].to(pos.unit).value / self._header["num_particles_total_side"]
		
		row = np.array([1,self._header["num_particles_total_side"],self._header["num_particles_total_side"]**2])
		posID = np.dot(pos.value/grid_unit,row)

		return posID 


	def visualize(self,fig=None,ax=None,scale=False,first=None,last=None,**kwargs):

		"""
		Visualize the particles in the Gadget snapshot using the matplotlib 3D plotting engine, the kwargs are passed to the matplotlib scatter method

		:param scale: if True, multiply all the (comoving) positions by the scale factor
		:type scale: bool.

		"""

		if not matplotlib:
			raise ImportError("matplotlib is not installed, cannot visualize!")

		#Get the positions if you didn't do it before
		if not hasattr(self,"positions"):
			self.getPositions()

		#If first or last are not specified, show all the particles
		if first is None:
			first = 0

		if last is None:
			last = self.positions.shape[0]

		#Instantiate figure
		if (fig is None) or (ax is None):
			
			self.fig = plt.figure()
			self.ax = self.fig.add_subplot(111,projection="3d")

		else:

			self.fig = fig
			self.ax = ax

		#Put the particles in the figure
		if scale:
			self.ax.scatter(*(self.positions[first:last].transpose()*self._header["scale_factor"]),**kwargs)
		else:
			self.ax.scatter(*self.positions[first:last].transpose(),**kwargs)

		#Put the labels on the axes
		self.ax.set_xlabel(r"$x({0})$".format(self.positions.unit.to_string()))
		self.ax.set_ylabel(r"$y({0})$".format(self.positions.unit.to_string()))
		self.ax.set_zlabel(r"$z({0})$".format(self.positions.unit.to_string()))

	def savefig(self,filename):

		"""
		Save the snapshot visualization to an external file

		:param filename: file name to which the figure will be saved
		:type filename: str.

		"""

		self.fig.savefig(filename)

	def close(self):

		"""
		Closes the snapshot file

		"""

		self.fp.close()

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
		_header_bare["masses"] = _header_bare["masses"].to(g).value * _header_bare["h"] / self._mass_unit
		_header_bare["num_particles_file_of_type"] = _header_bare["num_particles_file_of_type"].astype(np.int32)
		_header_bare["num_particles_total_of_type"] = _header_bare["num_particles_total_of_type"].astype(np.int32)
		_header_bare["comoving_distance"] = _header_bare["comoving_distance"].to(self.Mpc_over_h).value * 1.0e3


		#Convert units for positions and velocities
		_positions_converted = self.positions.to(self.kpc_over_h).value.astype(np.float32)

		if hasattr(self,"velocities"):
			assert self.positions.shape==self.velocities.shape
			_velocities_converted = (self.velocities.to(cm/s).value / self._velocity_unit).astype(np.float32)
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
				ext._gadget.write(_header_bare,_positions_converted[n*particles_per_file:(n+1)*particles_per_file],_velocities_converted[n*particles_per_file:(n+1)*particles_per_file],n*particles_per_file+1,filename_with_extension,writeVel)

			
			#The last file might have a different number of particles
			particles_last_file = len(_positions_converted[particles_per_file*(files-1):])

			#Update header
			_header_bare["num_particles_file"] = particles_last_file
			#TODO all particles are DM, fix distribution in the future
			_header_bare["num_particles_file_of_type"] = np.array([0,particles_last_file,0,0,0,0],dtype=np.int32)

			#Write it!
			filename_with_extension = "{0}.{1}".format(filename,files-1)
			ext._gadget.write(_header_bare,_positions_converted[particles_per_file*(files-1):],_velocities_converted[particles_per_file*(files-1):],(files-1)*particles_per_file+1,filename_with_extension,writeVel)

		else:

			#Update the num_files key only if not present already
			if "num_files" not in _header_bare.keys():
				_header_bare["num_files"] = 1

			#Update the header with the file names
			self.header["files"] = [ filename ]
			
			#Write it!!
			ext._gadget.write(_header_bare,_positions_converted,_velocities_converted,1,filename,writeVel)


	def writeParameterFile(self,filename,settings=Gadget2Settings.default()):

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


	
	def setPositions(self,positions):

		"""
		Sets the positions in the current snapshot (with the intent of writing them to a properly formatted snapshot file)

		:param positions: positions of the particles, must have units
		:type positions: (N,3) array with units

		"""

		assert positions.shape[1]==3
		assert positions.unit.physical_type=="length"

		self.positions = positions

	def setVelocities(self,velocities):

		"""
		Sets the velocities in the current snapshot (with the intent of writing them to a properly formatted snapshot file)

		:param velocities: velocities of the particles, must have units
		:type velocities: (N,3) array with units

		"""

		assert velocities.shape[1]==3
		assert velocities.unit.physical_type=="speed"

		self.velocities = velocities

	
	def setHeaderInfo(self,Om0=0.26,Ode0=0.74,w0=-1.0,wa=0.0,h=0.72,redshift=100.0,box_size=15.0*Mpc/0.72,flag_cooling=0,flag_sfr=0,flag_feedback=0,flag_stellarage=0,flag_metals=0,flag_entropy_instead_u=0,masses=np.array([0,1.03e10,0,0,0,0])*Msun,num_particles_file_of_type=None,npartTotalHighWord=np.zeros(6,dtype=np.uint32)):

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
		self.kpc_over_h = def_unit("kpc/h",kpc/self._header["h"])
		self.Mpc_over_h = def_unit("Mpc/h",Mpc/self._header["h"])

		#Compute the comoving distance according to the model
		cosmo = w0waCDM(H0=100.0*h,Om0=Om0,Ode0=Ode0,w0=w0,wa=wa)
		self._header["comoving_distance"] = cosmo.comoving_distance(redshift).to(self.Mpc_over_h)


	def numberDensity(self,resolution=0.5*Mpc,smooth=None,left_corner=None,save=False):

		"""
		Uses a C backend gridding function to compute the particle number density for the current snapshot: the density is evaluated using a nearest neighbor search

		:param resolution: resolution below which particles are grouped together; if an int is passed, this is the size of the grid
		:type resolution: float with units or int.

		:param smooth: if not None, performs a smoothing of the density (or potential) with a gaussian kernel of scale "smooth x the pixel resolution"
		:type smooth: int. or None

		:param left_corner: specify the position of the lower left corner of the box; if None, the minimum of the (x,y,z) of the contained particles is assumed
		:type left_corner: tuple of quantities or None

		:param save: if True saves the density histogram and resolution as instance attributes
		:type save: bool.

		:returns: tuple(numpy 3D array with the (unsmoothed) particle number density on a grid,bin resolution along the axes)  

		"""

		#Sanity checks
		assert type(resolution) in [np.int,quantity.Quantity]
		
		if type(resolution)==quantity.Quantity:	
			assert resolution.unit.physical_type=="length"

		#Check if positions are already available, otherwise retrieve them
		if hasattr(self,"positions"):
			positions = self.positions.copy()
		else:
			positions = self.getPositions(save=False)

		#Bin extremes (we start from the leftmost position up to the box size)
		if left_corner is None:
			xmin,ymin,zmin = positions.min(axis=0)
		else:
			xmin,ymin,zmin = left_corner

		#Construct binning
		if type(resolution)==quantity.Quantity:

			#Scale to appropriate units
			resolution = resolution.to(positions.unit)
			xi = np.arange(xmin.to(positions.unit).value,(xmin + self._header["box_size"]).to(positions.unit).value,resolution.value)
			yi = np.arange(ymin.to(positions.unit).value,(ymin + self._header["box_size"]).to(positions.unit).value,resolution.value)
			zi = np.arange(zmin.to(positions.unit).value,(zmin + self._header["box_size"]).to(positions.unit).value,resolution.value)

		else:

			xi = np.linspace(xmin.to(positions.unit).value,(xmin + self._header["box_size"]).to(positions.unit).value,resolution+1)
			yi = np.linspace(ymin.to(positions.unit).value,(ymin + self._header["box_size"]).to(positions.unit).value,resolution+1)
			zi = np.linspace(zmin.to(positions.unit).value,(zmin + self._header["box_size"]).to(positions.unit).value,resolution+1)


		#Compute the number count histogram
		assert positions.value.dtype==np.float32
		density = ext._gadget.grid3d(positions.value,(xi,yi,zi))

		#Accumulate from the other processors
		if self.pool is not None:
			
			self.pool.openWindow(density)
			self.pool.accumulate()
			self.pool.closeWindow()

		#Recompute resolution to make sure it represents the bin size correctly
		bin_resolution = ((xi[1:]-xi[:-1]).mean() * positions.unit,(yi[1:]-yi[:-1]).mean() * positions.unit,(zi[1:]-zi[:-1]).mean() * positions.unit)

		#Perform smoothing if prompted
		if smooth is not None:

			#Fourier transform the density field
			fx,fy,fz = np.meshgrid(fftfreq(density.shape[0]),fftfreq(density.shape[1]),rfftfreq(density.shape[2]),indexing="ij")
			density_ft = rfftn(density)

			#Perform the smoothing
			density_ft *= np.exp(-0.5*((2.0*np.pi*smooth)**2)*(fx**2 + fy**2 + fz**2))

			#Go back in real space
			density = irfftn(density_ft)


		#Return the density histogram, along with the bin resolution along each axis
		if save:
			self.density,self.resolution = density,bin_resolution

		return density,bin_resolution

	###################################################################################################################################################

	def cutPlaneGaussianGrid(self,normal=2,thickness=0.5*Mpc,center=7.0*Mpc,plane_resolution=0.1*Mpc,left_corner=None,thickness_resolution=0.1*Mpc,smooth=None,tomography=False,kind="density"):

		"""
		Cuts a density (or gravitational potential) plane out of the snapshot by computing the particle number density on a slab and performing Gaussian smoothing; the plane coordinates are cartesian comoving

		:param normal: direction of the normal to the plane (0 is x, 1 is y and 2 is z)
		:type normal: int. (0,1,2)

		:param thickness: thickness of the plane
		:type thickness: float. with units

		:param center: location of the plane along the normal direction
		:type center: float. with units

		:param plane_resolution: plane resolution (perpendicular to the normal)
		:type plane_resolution: float. with units (or int.)

		:param left_corner: specify the position of the lower left corner of the box; if None, the minimum of the (x,y,z) of the contained particles is assumed
		:type left_corner: tuple of quantities or None

		:param thickness_resolution: plane resolution (along the normal)
		:type thickness_resolution: float. with units (or int.)

		:param smooth: if not None, performs a smoothing of the density (or potential) with a gaussian kernel of scale "smooth x the pixel resolution"
		:type smooth: int. or None

		:param kind: decide if computing a density or gravitational potential plane (this is computed solving the poisson equation)
		:type kind: str. ("density" or "potential")

		:returns: tuple(numpy 3D array with the (unsmoothed) particle number density,bin resolution along the axes, number of particles on the plane)

		"""

		#Sanity checks
		assert normal in range(3),"There are only 3 dimensions!"
		assert kind in ["density","potential"],"Specify density or potential plane!"
		assert type(thickness)==quantity.Quantity and thickness.unit.physical_type=="length"
		assert type(center)==quantity.Quantity and center.unit.physical_type=="length"

		#Cosmological normalization factor
		cosmo_normalization = 1.5 * self._header["H0"]**2 * self._header["Om0"] / c**2

		#Direction of the plane
		plane_directions = range(3)
		plane_directions.pop(normal)

		#Get the particle positions if not available get
		if hasattr(self,"positions"):
			positions = self.positions.copy()
		else:
			positions = self.getPositions(save=False)

		#Lower left corner of the plane
		if left_corner is None:
			left_corner = positions.min(axis=0)

		#Create a list that holds the bins
		binning = [None,None,None]
		
		#Binning in the longitudinal direction
		assert type(plane_resolution) in [np.int,quantity.Quantity]
		
		if type(plane_resolution)==quantity.Quantity:
			
			assert plane_resolution.unit.physical_type=="length"
			plane_resolution = plane_resolution.to(positions.unit)
			binning[plane_directions[0]] = np.arange(left_corner[plane_directions[0]].to(positions.unit).value,(left_corner[plane_directions[0]] + self._header["box_size"]).to(positions.unit).value,plane_resolution.value)
			binning[plane_directions[1]] = np.arange(left_corner[plane_directions[1]].to(positions.unit).value,(left_corner[plane_directions[1]] + self._header["box_size"]).to(positions.unit).value,plane_resolution.value)

		else:

			binning[plane_directions[0]] = np.linspace(left_corner[plane_directions[0]].to(positions.unit).value,(left_corner[plane_directions[0]] + self._header["box_size"]).to(positions.unit).value,plane_resolution+1)
			binning[plane_directions[1]] = np.linspace(left_corner[plane_directions[1]].to(positions.unit).value,(left_corner[plane_directions[1]] + self._header["box_size"]).to(positions.unit).value,plane_resolution+1)

		
		#Binning in the normal direction		
		assert type(thickness_resolution) in [np.int,quantity.Quantity]
		center = center.to(positions.unit)
		thickness  = thickness.to(positions.unit)
		
		if type(thickness_resolution)==quantity.Quantity:
			
			assert thickness_resolution.unit.physical_type=="length"
			thickness_resolution = thickness_resolution.to(positions.unit)
			binning[normal] = np.arange((center - thickness/2).to(positions.unit).value,(center + thickness/2).to(positions.unit).value,thickness_resolution.value)

		else:

			binning[normal] = np.linspace((center - thickness/2).to(positions.unit).value,(center + thickness/2).to(positions.unit).value,thickness_resolution+1)

		#Now use gridding to compute the density along the slab
		assert positions.value.dtype==np.float32
		density = ext._gadget.grid3d(positions.value,tuple(binning))

		#Accumulate the density from the other processors
		if self.pool is not None:
			
			self.pool.openWindow(density)
			self.pool.accumulate()
			self.pool.closeWindow()

		#Compute the number of particles on the plane
		NumPartTotal = density.sum()

		#Recompute resolution to make sure it represents the bin size correctly
		bin_resolution = [(binning[0][1:]-binning[0][:-1]).mean() * positions.unit,(binning[1][1:]-binning[1][:-1]).mean() * positions.unit,(binning[2][1:]-binning[2][:-1]).mean() * positions.unit]

		#Longitudinal normalization factor
		density_normalization = bin_resolution[normal] * self._header["comoving_distance"] / self._header["scale_factor"]

		#Normalize the density to the density fluctuation
		density /= self._header["num_particles_total"]
		density *= (self._header["box_size"]**3 / (bin_resolution[0]*bin_resolution[1]*bin_resolution[2])).decompose().value

		#########################################################################################################################################
		######################################Decide if returning full tomographic information###################################################
		#########################################################################################################################################

		if tomography:

			if kind=="potential":
				raise NotImplementedError("Lensing potential tomography is not implemented!")

			if smooth is not None:
		
				#Fourier transform the density field
				fx,fy,fz = np.meshgrid(fftfreq(density.shape[0]),fftfreq(density.shape[1]),rfftfreq(density.shape[2]),indexing="ij")
				density_ft = rfftn(density)

				#Perform the smoothing
				density_ft *= np.exp(-0.5*((2.0*np.pi*smooth)**2)*(fx**2 + fy**2 + fz**2))

				#Go back in real space
				density = irfftn(density_ft)

			#Return the computed density histogram
			return density,bin_resolution,NumPartTotal


		#################################################################################################################################
		######################################Ready to solve poisson equation via FFTs###################################################
		#################################################################################################################################

		#First project along the normal direction
		density = density.sum(normal)
		bin_resolution.pop(normal)

		#If smoothing is enabled or potential calculations are needed, we need to FFT the density field
		if (smooth is not None) or kind=="potential":

			#Compute the multipoles
			lx,ly = np.meshgrid(fftfreq(density.shape[0]),rfftfreq(density.shape[1]),indexing="ij")
			l_squared = lx**2 + ly**2

			#Avoid dividing by 0
			l_squared[0,0] = 1.0

			#FFT the density field
			density_ft = rfftn(density)

			#Zero out the zeroth frequency
			density_ft[0,0] = 0.0

			if kind=="potential":
				#Solve the poisson equation
				density_ft *= -2.0 * (bin_resolution[0] * bin_resolution[1] / self.header["comoving_distance"]**2).decompose().value / (l_squared * ((2.0*np.pi)**2))

			if smooth is not None:
				#Perform the smoothing
				density_ft *= np.exp(-0.5*((2.0*np.pi*smooth)**2)*l_squared)

			#Revert the FFT
			density = irfftn(density_ft)

		#Multiply by the normalization factors
		density = density * cosmo_normalization * density_normalization
		density = density.decompose()
		assert density.unit.physical_type=="dimensionless"

		if kind=="potential":
			density *= rad**2
		else:
			density = density.value

		#Return
		return density,bin_resolution,NumPartTotal


	############################################################################################################################################################################

	def neighborDistances(self,neighbors=64):

		"""
		Find the N-th nearest neighbors to each particle

		:param neighbors: neighbor order
		:type neighbors: int.

		:returns: array with units

		"""

		#Get the particle positions if not available get
		if hasattr(self,"positions"):
			positions = self.positions.copy()
		else:
			positions = self.getPositions(save=False)

		#Build the KD-Tree
		particle_tree = KDTree(positions.value)

		#For memory reasons, with large datasets it's better to proceed in chunks with nearest neighbors queries
		numPart = positions.shape[0]
		rp = np.zeros(numPart)

		#Split the particles in chunks
		chunkSize = numPart // neighbors
		remaining = numPart % neighbors

		#Cycle over the chunks, querying the tree
		for i in range(neighbors):
			rp[i*chunkSize:(i+1)*chunkSize] = particle_tree.query(positions[i*chunkSize:(i+1)*chunkSize].value,k=neighbors)[0][:,neighbors-1]

		if remaining:
			rp[neighbors*chunkSize:] = particle_tree.query(positions[neighbors*chunkSize:].value,k=neighbors)[0][:,neighbors-1]

		#Return
		return rp * positions.unit


	############################################################################################################################################################################

	def cutPlaneAdaptive(self,normal=2,center=7.0*Mpc,left_corner=None,plane_resolution=0.1*Mpc,neighbors=64,neighborDistances=None,kind="density",projectAll=False):

		"""
		Cuts a density (or gravitational potential) plane out of the snapshot by computing the particle number density using an adaptive smoothing scheme; the plane coordinates are cartesian comoving

		:param normal: direction of the normal to the plane (0 is x, 1 is y and 2 is z)
		:type normal: int. (0,1,2)

		:param center: location of the plane along the normal direction
		:type center: float. with units

		:param plane_resolution: plane resolution (perpendicular to the normal)
		:type plane_resolution: float. with units (or int.)

		:param left_corner: specify the position of the lower left corner of the box; if None, the minimum of the (x,y,z) of the contained particles is assumed
		:type left_corner: tuple of quantities or None

		:param neighbors: number of nearest neighbors to use in the adaptive smoothing procedure
		:type neighbors: int.

		:param neighborDistances: precomputed distances of each particle to its N-th nearest neighbor; if None these are computed
		:type neighborDistances: array with units

		:param kind: decide if computing a density or gravitational potential plane (this is computed solving the poisson equation)
		:type kind: str. ("density" or "potential")

		:param projectAll: if True, all the snapshot is projected on a single slab perpendicular to the normal, ignoring the position of the center
		:type projectAll: bool.

		:returns: tuple(numpy 2D array with the computed particle number density (or lensing potential),bin resolution along the axes,number of particles on the plane)

		"""

		#Sanity checks
		assert normal in range(3),"There are only 3 dimensions!"
		assert kind in ["density","potential"],"Specify density or potential plane!"
		assert type(center)==quantity.Quantity and center.unit.physical_type=="length"

		#Direction of the plane
		plane_directions = range(3)
		plane_directions.pop(normal)

		#Get the particle positions if not available get
		if hasattr(self,"positions"):
			positions = self.positions.copy()
		else:
			positions = self.getPositions(save=False)

		#Lower left corner of the plane
		if left_corner is None:
			left_corner = positions.min(axis=0)
		
		#Binning of the plane
		binning = [None,None]
		assert type(plane_resolution) in [np.int,quantity.Quantity]
		
		if type(plane_resolution)==quantity.Quantity:
			
			assert plane_resolution.unit.physical_type=="length"
			plane_resolution = plane_resolution.to(positions.unit)

			for i in range(2):
				binning[i] = np.arange(left_corner[plane_directions[i]].to(positions.unit).value,(left_corner[plane_directions[i]] + self._header["box_size"]).to(positions.unit).value,plane_resolution.value)

		else:

			for i in range(2):
				binning[i] = np.linspace(left_corner[plane_directions[i]].to(positions.unit).value,(left_corner[plane_directions[i]] + self._header["box_size"]).to(positions.unit).value,plane_resolution+1)

		#Recompute bin_resolution
		bin_resolution = [ (binning[0][1:]-binning[0][:-1]).mean() * positions.unit,(binning[1][1:]-binning[1][:-1]).mean() * positions.unit ]

		###################################################################################
		#For each particle, we need to determine the distance to its N-th nearest neighbor#
		###################################################################################

		if neighborDistances is None:
	
			#Find the distance to the Nth-nearest neighbor
			rp = self.neighborDistances(neighbors).to(positions.unit).value

		else:
			
			#Convert pre computed distances into appropriate units
			assert neighbors is None,"You cannot specify the number of neighbors if the distances are precomputed!"
			assert neighborDistances.shape[0]==positions.shape[0]
			rp = neighborDistances.to(positions.unit).value

		#Check that thay are all positive
		assert (rp>0).all()

		#Compute the adaptive smoothing
		density = (3.0/np.pi)*ext._gadget.adaptive(positions.value,rp,binning,center.to(positions.unit).value,plane_directions[0],plane_directions[1],normal,projectAll)

		#Accumulate the density from the other processors
		if self.pool is not None:
			
			self.pool.openWindow(density)
			self.pool.accumulate()
			self.pool.closeWindow()

		#Integrate the density to find the total number of particles
		NumPartTotal = (density.sum() * bin_resolution[0] * bin_resolution[1] * positions.unit**-2).decompose().value

		##############################################
		#Compute the dimensionless density fluctation#
		##############################################

		#Normalize to correct units and subtract the mean
		density *= positions.unit**-2
		density *= (self.header["box_size"]**3 / self.header["num_particles_total"]).decompose()
		density -= self.header["box_size"]

		#Add the cosmological normalization factor
		density *= 1.5 * self.header["H0"]**2 * self.header["Om0"] / c**2
		density *= self.header["comoving_distance"] / self.header["scale_factor"]
		assert density.unit.physical_type=="dimensionless" 
		density = density.decompose().value

		if kind=="density":
			return density,bin_resolution,NumPartTotal

		#################################################################################
		##############Ready to compute the lensing potential#############################
		#################################################################################

		if kind=="potential":

			#Compute the multipoles
			lx,ly = np.meshgrid(fftfreq(density.shape[0]),rfftfreq(density.shape[1]),indexing="ij")
			l_squared = lx**2 + ly**2

			#Avoid dividing by 0
			l_squared[0,0] = 1.0

			#FFT the density field
			density_ft = rfftn(density)

			#Zero out the zeroth frequency
			density_ft[0,0] = 0.0

			#Solve the poisson equation
			density_ft *= -2.0 * (bin_resolution[0] * bin_resolution[1] / self.header["comoving_distance"]**2).decompose().value / (l_squared * ((2.0*np.pi)**2))

			#Revert the FFT and return
			density = irfftn(density_ft)
			return density*(rad**2),bin_resolution,NumPartTotal



	############################################################################################################################################################################

	def cutPlaneAngular(self,normal=2,thickness=0.5*Mpc,center=7.0*Mpc,left_corner=None,plane_lower_corner=np.array([0.0,0.0])*deg,plane_size=0.15*deg,plane_resolution=1.0*arcmin,thickness_resolution=0.1*Mpc,smooth=None,tomography=False,kind="density",space="real"):

		"""
		Same as cutPlaneGaussianGrid(), except that this method will return a lens plane as seen from an observer at z=0; the spatial transverse units are converted in angular units as seen from the observer

		:param normal: direction of the normal to the plane (0 is x, 1 is y and 2 is z)
		:type normal: int. (0,1,2)

		:param thickness: thickness of the plane
		:type thickness: float. with units

		:param center: location of the plane along the normal direction; it is assumed that the center of the plane is seen from an observer with a redshift of self.header["redshift"]
		:type center: float. with units

		:param left_corner: specify the position of the lower left corner of the box; if None, the minimum of the (x,y,z) of the contained particles is assumed
		:type left_corner: tuple of quantities or None

		:param plane_lower_corner: lower left corner of the plane, as seen from the observer (0,0) corresponds to the lower left corner of the snapshot
		:type plane_lower_corner: float with units.

		:param plane_size: angular size of the lens plane (angles start from 0 in the lower left corner)
		:type plane_size: float with units

		:param plane_resolution: plane angular resolution (perpendicular to the normal)
		:type plane_resolution: float. with units (or int.)

		:param thickness_resolution: plane resolution (along the normal)
		:type thickness_resolution: float. with units (or int.)

		:param smooth: if not None, performs a smoothing of the angular density (or potential) with a gaussian kernel of scale "smooth x the pixel resolution"
		:type smooth: int. or None

		:param tomography: if True returns the lens plane angular density for each slab, otherwise a projected density (or lensing potential) is computed
		:type tomography: bool.

		:param kind: decide if computing an angular density or lensing potential plane (this is computed solving the poisson equation)
		:type kind: str. ("density" or "potential")

		:param space: if "real" return the lens plane in real space, if "fourier" the Fourier transform is not inverted
		:type space: str.

		:returns: tuple(numpy 2D or 3D array with the (unsmoothed) particle angular number density,bin angular resolution, total number of particles on the plane); the constant spatial part of the density field is subtracted (we keep the fluctuation only)

		"""

		#Sanity checks
		assert normal in range(3),"There are only 3 dimensions!"
		assert kind in ["density","potential"],"Specify density or potential plane!"
		assert type(thickness)==quantity.Quantity and thickness.unit.physical_type=="length"
		assert type(center)==quantity.Quantity and center.unit.physical_type=="length"
		assert type(plane_lower_corner)==quantity.Quantity and plane_lower_corner.unit.physical_type=="angle"
		assert type(plane_size)==quantity.Quantity and plane_size.unit.physical_type=="angle"

		#First compute the overall normalization factor for the angular density
		cosmo_normalization = 1.5 * (self._header["H0"]**2) * self._header["Om0"]  * self.cosmology.comoving_distance(self._header["redshift"]) * (1.0+self._header["redshift"]) / c**2

		#Direction of the plane
		plane_directions = range(3)
		plane_directions.pop(normal)

		#Get the particle positions if not available get
		if hasattr(self,"positions"):
			positions = self.positions.copy()
		else:
			positions = self.getPositions(save=False)

		#Scale the units
		thickness = thickness.to(positions.unit)
		center = center.to(positions.unit)

		#Lower left corner of the plane
		if left_corner is None:
			left_corner = positions.min(axis=0)

		#Translate the transverse coordinates so that the lower corner is in (0,0)
		for i in range(2):
			positions[:,plane_directions[i]] -= left_corner[plane_directions[i]].astype(np.float32)

		#Create a list that holds the bins
		binning = [None,None,None]
		
		#Binning in the longitudinal direction
		assert type(plane_resolution) in [np.int,quantity.Quantity]
		
		if type(plane_resolution)==quantity.Quantity:
			
			assert plane_resolution.unit.physical_type=="angle"
			plane_resolution = plane_resolution.to(rad)
			binning[plane_directions[0]] = np.arange(plane_lower_corner[0].to(rad).value,(plane_lower_corner[0] + plane_size).to(rad).value,plane_resolution.value)
			binning[plane_directions[1]] = np.arange(plane_lower_corner[1].to(rad).value,(plane_lower_corner[1] + plane_size).to(rad).value,plane_resolution.value)

		else:

			binning[plane_directions[0]] = np.linspace(plane_lower_corner[0].to(rad).value,(plane_lower_corner[0] + plane_size).to(rad).value,plane_resolution + 1)
			binning[plane_directions[1]] = np.linspace(plane_lower_corner[1].to(rad).value,(plane_lower_corner[1] + plane_size).to(rad).value,plane_resolution + 1)

		
		#Get the snapshot comoving distance from the observer (which is the same as the plane comoving distance)
		plane_comoving_distance = self.cosmology.comoving_distance(self._header["redshift"]).to(positions.unit)

		#Binning in the normal direction		
		assert type(thickness_resolution) in [np.int,quantity.Quantity]
		center = center.to(positions.unit)
		thickness  = thickness.to(positions.unit)
		
		if type(thickness_resolution)==quantity.Quantity:
			
			assert thickness_resolution.unit.physical_type=="length"
			thickness_resolution = thickness_resolution.to(positions.unit)
			binning[normal] = np.arange((plane_comoving_distance - thickness/2).to(positions.unit).value,(plane_comoving_distance + thickness/2).to(positions.unit).value,thickness_resolution.value)

		else:

			binning[normal] = np.linspace((plane_comoving_distance - thickness/2).to(positions.unit).value,(plane_comoving_distance + thickness/2).to(positions.unit).value,thickness_resolution+1)


		#Now that everything has the same units, let's go dimensionless to convert into angular units
		length_unit = positions.unit
		positions = positions.value

		#Convert the normal direction into comoving distance from the observer
		positions[:,normal] += (plane_comoving_distance.value - center.value)

		#Convert the longitudinal spatial coordinates into angles (theta = comiving transverse/comoving distance)
		for i in range(2):
			positions[:,plane_directions[i]] /= positions[:,normal]

		#Now use grid3d to compute the angular density on the lens plane
		assert positions.dtype==np.float32
		density = ext._gadget.grid3d(positions,tuple(binning))

		#Accumulate the density from the other processors
		if self.pool is not None:
			
			self.pool.openWindow(density)
			self.pool.accumulate()
			self.pool.closeWindow()


		#Compute the total number of particles on the lens plane
		NumPartTotal = density.sum()

		#Recompute resolution to make sure it represents the bin size correctly
		bin_resolution = [ (binning[0][1:]-binning[0][:-1]).mean() , (binning[1][1:]-binning[1][:-1]).mean() , (binning[2][1:]-binning[2][:-1]).mean() ]
		
		#Restore units
		bin_resolution[normal] *= length_unit
	
		for i in range(2):

			try:
				bin_resolution[plane_directions[i]] = (bin_resolution[plane_directions[i]] * rad).to(plane_resolution.unit)
			except AttributeError:
				bin_resolution[plane_directions[i]] = (bin_resolution[plane_directions[i]] * rad).to(arcmin) 

		#############################################################################################################################################
		######################################If tomography is desired, we can return now############################################################
		#############################################################################################################################################

		if tomography:

			if kind=="potential":
				raise NotImplementedError("Lensing potential tomography is not implemented!")

			if smooth is not None:

				fx,fy,fz = np.meshgrid(fftfreq(density.shape[0]),fftfreq(density.shape[1]),rfftfreq(density.shape[2]),indexing="ij")
				density_ft = rfftn(density)
				density_ft *= np.exp(-0.5*((2.0*np.pi*smooth)**2)*(fx**2 + fy**2 + fz**2))
				density_ft[0,0] = 0.0
				density = irfftn(density_ft)

				return (density * (1.0/self._header["num_particles_total"]) * (self._header["box_size"]*self.lensMaxSize()**2)/reduce(mul,bin_resolution)).decompose().value, bin_resolution, NumPartTotal

			else:

				return ((density - density.sum()/reduce(mul,density.shape)) * (1.0/self._header["num_particles_total"]) * (self._header["box_size"]*self.lensMaxSize()**2)/reduce(mul,bin_resolution)).decompose().value, bin_resolution, NumPartTotal

		#############################################################################################################################################
		######################################Ready to solve the lensing poisson equation via FFTs###################################################
		#############################################################################################################################################

		#First project the density along the line of sight
		density = density.sum(normal)
		bin_resolution.pop(normal)

		#Compute the normalization factor to convert the absolute number density into a relative number density
		density_normalization = (self._header["box_size"]/self._header["num_particles_total"]) * (self.lensMaxSize() / bin_resolution[0])**2

		#Then solve the poisson equation and/or smooth the density field with FFTs
		if (smooth is not None) or kind=="potential":
		
			#Compute the multipoles
			lx,ly = np.meshgrid(fftfreq(density.shape[0]),rfftfreq(density.shape[1]),indexing="ij")
			l_squared = lx**2 + ly**2

			#Avoid dividing by 0
			l_squared[0,0] = 1.0

			#Fourier transform the density field
			density_ft = rfftn(density)

			#Perform the smoothing
			if smooth is not None:
				density_ft *= np.exp(-0.5*((2.0*np.pi*smooth)**2)*l_squared)

			#If kind is potential, solve the poisson equation
			if kind=="potential":
				density_ft *= -2.0 * ((bin_resolution[0].to(rad).value)**2) / (l_squared * ((2.0*np.pi)**2))
				
			#Return only the density fluctuation, dropping the zeroth frequency (i.e. uniform part)
			density_ft[0,0] = 0.0 

			#Go back in real space
			if space=="real":
				density = irfftn(density_ft)
			elif space=="fourier":
				density = density_ft
			else:
				raise ValueError("space must be real or fourier!")

		else:

			density -= density.sum() / reduce(mul,density.shape)
			if space=="fourier":
				density = rfftn(density)

		#Return
		return (density*cosmo_normalization*density_normalization).decompose().value,bin_resolution,NumPartTotal


	#############################################################################################################################################

	def lensMaxSize(self):

		"""
		Computes the maximum observed size of a lens plane cut out of the current snapshot
	
		"""

		return ((self._header["box_size"] / self.cosmology.comoving_distance(self._header["redshift"])) * rad).to(deg)


	#############################################################################################################################################


	def powerSpectrum(self,k_edges,resolution=None):

		"""
		Computes the power spectrum of the relative density fluctuations in the snapshot at the wavenumbers specified by k_edges; a discrete particle number density is computed before hand to prepare the FFT grid

		:param k_edges: wavenumbers at which to compute the density power spectrum (must have units)
		:type k_edges: array.

		:param resolution: optional, fix the grid resolution to some value; to be passed to the numberDensity method. If none this is computed automatically from the k_edges
		:type resolution: float with units, int. or None

		:returns: tuple(k_values(bin centers),power spectrum at the specified k_values)

		"""

		#Check for correct units
		assert k_edges.unit.physical_type=="wavenumber"

		if resolution is None:
			resolution = 2.0 * np.pi / k_edges.max()

		#Sanity check on bin spacing (must not be smaller than the one allowed by the size of the box)
		if (k_edges[1:] - k_edges[:-1]).mean() < 2.0*np.pi/self._header["box_size"]:
			raise ValueError("Your bins are too small! Minimum allowed by the current box size is {0}".format(2.0*np.pi/self._header["box_size"]))

		#Compute the gridded number density
		if not hasattr(self,"density"):
			density,bin_resolution = self.numberDensity(resolution=resolution) 
		else:
			assert resolution is None,"The spatial resolution is already specified in the attributes of this instance! Call numberDensity() to modify!"
			density,bin_resolution = self.density,self.resolution
		
		#Decide pixel sizes in Fourier spaces
		kpixX = (2.0*np.pi/self._header["box_size"]).to(k_edges.unit)
		kpixY = (2.0*np.pi/self._header["box_size"]).to(k_edges.unit)
		kpixZ = (2.0*np.pi/self._header["box_size"]).to(k_edges.unit)

		#Compute the maximum allowed wavenumber
		k_max = 0.5*np.sqrt((kpixX * density.shape[0])**2 + (kpixY * density.shape[1])**2 + (kpixZ * density.shape[2])**2)
		k_max_recommended = (1 / (max(bin_resolution))).to(k_max.unit)

		#Sanity check on maximum k: maximum is limited by the grid resolution
		if k_edges.max() > k_max:
			print("Your grid resolution is too low to compute accurately the power on {0} (maximum recommended {1}, distortions might start to appear already at {2}): results might be inaccurate".format(k_edges.max(),k_max,k_max_recommended))

		#Perform the FFT
		density_ft = rfftn(density)

		#Compute the azimuthal averages
		hits,power_spectrum = ext._topology.rfft3_azimuthal(density_ft,density_ft,kpixX.value,kpixY.value,kpixZ.value,k_edges.value)

		#Return the result (normalize the power so it corresponds to the one of the density fluctuations)
		k = 0.5*(k_edges[1:]+k_edges[:-1])
		return k,(power_spectrum/hits) * (self._header["box_size"]**3) / (self._header["num_particles_total"]**2)



	def __add__(self,rhs):

		"""
		Add two gadget snapshots together: useful when the particle content is split between different files; all the positions and particle velocities are vstacked together

		"""

		merged_snapshot = Gadget2Snapshot(None)
		merged_snapshot._header = self._header + rhs._header

		if hasattr(self,"positions") and hasattr(rhs,"positions"):
			
			assert self.positions.unit==rhs.positions.unit
			merged_snapshot.positions = np.vstack((self.positions.value,rhs.positions.value))
			merged_snapshot.positions *= self.positions.unit

		if hasattr(self,"velocities") and hasattr(rhs,"velocities"):
			
			assert self.velocities.unit==rhs.velocities.unit
			merged_snapshot.velocities = np.vstack((self.velocities.value,rhs.velocities.value))
			merged_snapshot.velocities *= self.velocities.unit

		if hasattr(self,"id") and hasattr(rhs,"id"):

			merged_snapshot.id = np.hstack((self.id,rhs.id))


		return merged_snapshot



