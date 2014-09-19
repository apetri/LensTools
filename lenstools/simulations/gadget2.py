from __future__ import division

from .. import external as ext

import numpy as np
from scipy.stats import rankdata
from astropy.units import kpc,Mpc,cm,km,g,s,Msun,quantity,def_unit

#FFT engines
from numpy.fft import rfftn

try:
	
	import matplotlib.pyplot as plt
	from mpl_toolkits.mplot3d import Axes3D
	matplotlib = True

except ImportError:

	matplotlib = False

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
		fields_to_match = ["Ode0","Om0","h","box_size","endianness","flag_cooling","flag_feedback","flag_sfr","num_files"]
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
	A class that handles I/O from Gadget2 snapshots: it was mainly designed to parse physical information from binary gadget snapshots

	"""

	_mass_unit = 1.989e43
	_velocity_unit = 1.0e5

	def __init__(self,fp=None):

		assert (type(fp)==file) or (fp is None),"Call the open() method instead!!"

		if fp is not None:
		
			self.fp = fp
			self._header = Gadget2Header(ext._gadget.getHeader(fp))
			self._header["files"] = [self.fp.name]

			#Scale box to kpc
			self._header["box_size"] *= kpc
			#Convert to Mpc
			self._header["box_size"] = self._header["box_size"].to(Mpc)

			#Scale masses to correct units
			self._header["masses"] *= (self._mass_unit / self._header["h"])
			self._header["masses"] *= g
			self._header["masses"] = self._header["masses"].to(Msun) 

			#Scale Hubble parameter to correct units
			self._header["H0"] = self._header["h"] * 100 * km / (s*Mpc)

			#Update the dictionary with the number of particles per side
			self._header["num_particles_total_side"] = int(np.round(self._header["num_particles_total"]**(1/3)))

			#Define the Mpc/h unit for convenience
			self.Mpc_over_h = def_unit("Mpc/h",Mpc/self._header["h"])

	@classmethod
	def open(cls,filename):

		"""
		Opens a gadget snapshot at filename

		:param filename: file name of the gadget snapshot
		:type filename: str. or file.
		"""

		if type(filename)==str:
			fp = open(filename,"r")
		elif type(filename)==file:
			fp = filename
		else:
			raise TypeError("filename must be string or file!")
		
		return cls(fp)

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
		positions = (ext._gadget.getPosVel(self.fp,offset,numPart) * kpc).to(Mpc) 
		if save:
			self.positions = positions
		
		#Return
		return positions

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
		ranks = rankdata(self.id).astype(np.int32) - 1

		#Sort positions
		if hasattr(self,"positions"):
			
			assert self.positions.shape[0]==len(self.id)
			self.positions = self.positions[ranks]

		#Sort velocities
		if hasattr(self,"velocities"):

			assert self.velocities.shape[0]==len(self.id)
			self.velocities = self.velocities[ranks]

		#Finally sort IDs
		self.id.sort() 


	def visualize(self,fig=None,ax=None,scale=False,**kwargs):

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

		#Instantiate figure
		if (fig is None) or (ax is None):
			
			self.fig = plt.figure()
			self.ax = self.fig.add_subplot(111,projection="3d")

		else:

			self.fig = fig
			self.ax = ax

		#Put the particles in the figure
		if scale:
			self.ax.scatter(*(self.positions.transpose()*self._header["scale_factor"]),**kwargs)
		else:
			self.ax.scatter(*self.positions.transpose(),**kwargs)

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
		assert hasattr(self,"positions") and hasattr(self,"velocities"),"Positions and velocities must be specified!!"
		assert self.positions.shape[0]==self.velocities.shape[0]

		if not hasattr(self,"_header"):
			self.setHeaderInfo()	

		#Build a bare header based on the available info (need to convert units back to the Gadget ones)
		_header_bare = self._header.copy()
		_header_bare["box_size"] = _header_bare["box_size"].to(kpc).value
		_header_bare["masses"] = _header_bare["masses"].to(g).value * _header_bare["h"] / self._mass_unit
		_header_bare["num_particles_file_of_type"] = _header_bare["num_particles_file_of_type"].astype(np.int32)
		_header_bare["num_particles_total_of_type"] = _header_bare["num_particles_file_of_type"].astype(np.int32)

		#Convert units for positions and velocities
		_positions_converted = self.positions.to(kpc).value.astype(np.float32)
		_velocities_converted = (self.velocities.to(cm/s).value / self._velocity_unit).astype(np.float32)

		#Check if we want to split on multiple files (only DM particles supported so far for this feature)
		if files>1:

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
				ext._gadget.write(_header_bare,_positions_converted[n*particles_per_file:(n+1)*particles_per_file],_velocities_converted[n*particles_per_file:(n+1)*particles_per_file],n*particles_per_file+1,filename_with_extension)

			
			#The last file might have a different number of particles
			particles_last_file = len(_positions_converted[particles_per_file*(files-1):])

			#Update header
			_header_bare["num_particles_file"] = particles_last_file
			#TODO all particles are DM, fix distribution in the future
			_header_bare["num_particles_file_of_type"] = np.array([0,particles_last_file,0,0,0,0],dtype=np.int32)

			#Write it!
			filename_with_extension = "{0}.{1}".format(filename,files-1)
			ext._gadget.write(_header_bare,_positions_converted[particles_per_file*(files-1):],_velocities_converted[particles_per_file*(files-1):],(files-1)*particles_per_file+1,filename_with_extension)

		else:
			
			#Write it!!
			ext._gadget.write(_header_bare,_positions_converted,_velocities_converted,1,filename)

	
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

	
	def setHeaderInfo(self,Om0=0.26,Ode0=0.74,h=0.72,redshift=1.0,box_size=15.0*Mpc,flag_cooling=0,flag_sfr=0,flag_feedback=0,masses=np.array([0,1.03e10,0,0,0,0])*Msun,num_particles_file_of_type=None):

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
		self._header["h"] = h
		self._header["redshift"] = redshift
		self._header["scale_factor"] = 1.0 / (1.0 + redshift)
		self._header["box_size"] = box_size
		self._header["flag_cooling"] = flag_cooling
		self._header["flag_sfr"] = flag_sfr
		self._header["flag_feedback"] = flag_feedback
		self._header["masses"] = masses
		self._header["num_particles_file_of_type"] = num_particles_file_of_type
		self._header["num_particles_file"] = num_particles_file_of_type.sum()
		self._header["num_particles_total_of_type"] = num_particles_file_of_type
		self._header["num_particles_total"] = num_particles_file_of_type.sum()

		#Define the Mpc/h unit for convenience
		self.Mpc_over_h = def_unit("Mpc/h",Mpc/self._header["h"])


	def numberDensity(self,resolution=0.5*Mpc,save=False):

		"""
		Uses the np.histogramdd function to compute the particle number density for the current snapshot: the density is evaluated using a nearest neighbor search

		:param resolution: resolution below which particles are grouped together; if an int is passed, this is the size of the grid
		:type resolution: float with units or int.

		:param save: if True saves the density histogram and resolution as instance attributes
		:type save:

		:returns: tuple(numpy 3D array with the (unsmoothed) particle number density,bin resolution along the axes)  

		"""

		#Sanity checks
		assert type(resolution) in [np.int,quantity.Quantity]
		
		if type(resolution)==quantity.Quantity:	
			assert resolution.unit.physical_type=="length"

		#Check if positions are already available, otherwise retrieve them
		if hasattr(self,"positions"):
			positions = self.positions
		else:
			positions = self.getPositions(save=False)

		#Bin extremes (we start from the leftmost position up to the box size)
		xmin,ymin,zmin = positions.min(axis=0)

		#Construct binning
		if type(resolution)==quantity.Quantity:

			#Scale to appropriate units
			resolution = resolution.to(positions.unit)
			xi = np.arange(xmin.value,(xmin + self._header["box_size"]).value,resolution.value)
			yi = np.arange(ymin.value,(ymin + self._header["box_size"]).value,resolution.value)
			zi = np.arange(zmin.value,(zmin + self._header["box_size"]).value,resolution.value)

		else:

			xi = np.linspace(xmin.value,(xmin + self._header["box_size"]).value,resolution+1)
			yi = np.linspace(ymin.value,(ymin + self._header["box_size"]).value,resolution+1)
			zi = np.linspace(zmin.value,(zmin + self._header["box_size"]).value,resolution+1)


		#Compute the number count histogram
		density,bins = np.histogramdd(positions.value,(xi,yi,zi))

		#Recompute resolution to make sure it represents the bin size correctly
		bin_resolution = ((bins[0][1:]-bins[0][:-1]).mean() * positions.unit,(bins[1][1:]-bins[1][:-1]).mean() * positions.unit,(bins[2][1:]-bins[2][:-1]).mean() * positions.unit)

		#Return the density histogram, along with the bin resolution along each axis
		if save:
			self.density,self.resolution = density,bin_resolution

		return density,bin_resolution


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
		k_max_recommended = 1 / (max(bin_resolution))

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



