from __future__ import division

from abc import ABCMeta,abstractproperty,abstractmethod

from operator import mul
from functools import reduce

import sys,os

from .logs import logplanes,logstderr,peakMemory
import logging

from .. import extern as ext

import numpy as np

#astropy stuff, invaluable here
from astropy.units import Mbyte,kpc,Mpc,cm,km,g,s,hour,day,deg,arcmin,rad,Msun,quantity,def_unit
from astropy.constants import c
from astropy.cosmology import w0waCDM,z_at_value

#FFT engine
from ..utils.fft import NUMPYFFTPack
fftengine = NUMPYFFTPack()

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


###################################################################
#################NbodySnapshot abstract class######################
###################################################################

class NbodySnapshot(object):

	__metaclass__ = ABCMeta

	"""
	A class that handles Nbody simulation snapshots; it's an abstract class as I/O routines have to be specified

	"""

	#####################################################################
	######################Abstract methods###############################
	#####################################################################

	@abstractmethod
	def buildFilename(cls,root,pool,**kwargs):
		pass

	@abstractmethod
	def int2root(cls,name,n):
		pass	

	@abstractmethod
	def getHeader(self):
		pass

	@abstractmethod
	def setLimits(self):
		pass

	@abstractmethod
	def getPositions(self,first=None,last=None,save=True):
		pass

	@abstractmethod
	def getVelocities(self,first=None,last=None,save=True):
		pass

	@abstractmethod
	def getID(self,first=None,last=None,save=True):
		pass

	@abstractmethod
	def write(self,filename,files=1):
		pass

	###################################################################################
	######################Default, non--abstract methods###############################
	###################################################################################

	#Check that header has all required keys#
	_header_keys = ['redshift','scale_factor','comoving_distance','masses','num_particles_file','num_particles_total','box_size','num_files','Om0','Ode0','w0','wa','h']

	def _check_header(self):

		for key in self._header_keys:
			assert key in self._header,"Key {0} not loaded in header, please make sure that the getHeader method is configured to do that!".format(key)

	####################################################################################################################

	def __enter__(self):
		return self

	def __exit__(self,type,value,tb):
		self.fp.close()

	def __init__(self,fp=None,pool=None,length_unit=1.0*kpc,mass_unit=1.0e10*Msun,velocity_unit=1.0*km/s,header_kwargs=dict()):

		self.pool = pool

		self._length_unit = length_unit.to(cm).value
		self._mass_unit = mass_unit.to(g).value
		self._velocity_unit = velocity_unit.to(cm/s).value

		if fp is not None:
		
			self.fp = fp

			#Load the header
			self._header = self.getHeader(**header_kwargs)
			
			#Check that header has been loaded correctly
			self._check_header()

			#Hubble parameter
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
				if "comoving_distance" in self._header:
					self._header["comoving_distance"] = (self._header["comoving_distance"] / 1.0e3) * self.Mpc_over_h

			else:
				self._header["box_size"] *= kpc
				logging.debug("Warning! Hubble parameter h is zero!!")

			#Scale masses to correct units
			if h>0.0:
				self._header["masses"] *= (self._mass_unit / self._header["h"])
				self._header["masses"] = (self._header["masses"]*g).to(Msun) 

			#Scale Hubble parameter to correct units
			self._header["H0"] = self._header["h"] * 100 * km / (s*Mpc)

			#Update the dictionary with the number of particles per side
			self._header["num_particles_total_side"] = int(np.round(self._header["num_particles_total"]**(1/3)))

			#Once all the info is available, add a wCDM instance as attribute to facilitate the cosmological calculations
			if h>0.0:
				self.cosmology = w0waCDM(H0=self._header["H0"],Om0=self._header["Om0"],Ode0=self._header["Ode0"],w0=self._header["w0"],wa=self._header["wa"])

			#Set particle number limits that this instance will handle
			self.setLimits()

	@classmethod
	def open(cls,filename,pool=None,header_kwargs=dict(),**kwargs):

		"""
		Opens a snapshot at filename

		:param filename: file name of the snapshot
		:type filename: str. or file.

		:param pool: use to distribute the calculations on different processors; if not None, each processor takes care of one of the snapshot parts, appending as ".n" to the filename
		:type pool: MPIWhirlPool instance

		:param header_kwargs: keyword arguments to pass to the getHeader method
		:type header_kwargs: dict.

		:param kwargs: the keyword arguments are passed to buildFilename
		:type kwargs: dict.

		"""

		if hasattr(filename,"format"):

			fp = open(cls.buildFilename(filename,pool,**kwargs),"rb")
		
		elif hasattr(filename,"read"):
			
			if pool is not None:
				raise TypeError("Specifying file objects with MPIPools is not allowed!")
			fp = filename
		
		else:
			raise TypeError("filename type is {0}, must be string or file!".format(type(filename)))
		
		return cls(fp,pool,header_kwargs=header_kwargs)

	@property
	def header(self):

		"""
		Displays the snapshot header information

		:returns: the snapshot header information in dictionary form
		:rtype: dict.

		"""

		return self._header


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
		Visualize the particles in the snapshot using the matplotlib 3D plotting engine, the kwargs are passed to the matplotlib scatter method

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
			self.ax.scatter(*(self.positions[first:last].T.value*self._header["scale_factor"]),**kwargs)
		else:
			self.ax.scatter(*self.positions[first:last].T.value,**kwargs)

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


	def massDensity(self,resolution=0.5*Mpc,smooth=None,left_corner=None,save=False,density_placeholder=None):

		"""
		Uses a C backend gridding function to compute the matter mass density fluctutation for the current snapshot: the density is evaluated using a nearest neighbor search

		:param resolution: resolution below which particles are grouped together; if an int is passed, this is the size of the grid
		:type resolution: float with units or int.

		:param smooth: if not None, performs a smoothing of the density (or potential) with a gaussian kernel of scale "smooth x the pixel resolution"
		:type smooth: int. or None

		:param left_corner: specify the position of the lower left corner of the box; if None, the minimum of the (x,y,z) of the contained particles is assumed
		:type left_corner: tuple of quantities or None

		:param save: if True saves the density histogram and resolution as instance attributes
		:type save: bool.

		:param density placeholder: if not None, it is used as a fixed memory chunk for MPI communications of the density
		:type density_placeholder: array

		:returns: tuple(numpy 3D array with the (unsmoothed) matter density fluctuation on a grid,bin resolution along the axes)  

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

		assert hasattr(self,"weights")
		assert hasattr(self,"virial_radius")
		assert hasattr(self,"concentration")

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

		#Weights
		if self.weights is not None:
			weights = (self.weights * self._header["num_particles_total"] / ((len(xi) - 1) * (len(yi) - 1) * (len(zi) - 1))).astype(np.float32)
		else:
			weights = None

		if self.virial_radius is not None:
			rv = self.virial_radius.to(positions.unit).value
		else:
			rv = None

		density = ext._nbody.grid3d(positions.value,(xi,yi,zi),weights,rv,self.concentration) * (len(xi)-1) * (len(yi)-1) * (len(zi)-1) / self._header["num_particles_total"]

		#Accumulate from the other processors
		if self.pool is not None:

			if density_placeholder is not None:

				density_placeholder[:] = density
				self.pool.comm.Barrier()
				self.pool.accumulate()

			else:
			
				self.pool.openWindow(density)
				self.pool.accumulate()
				self.pool.closeWindow()

		#Recompute resolution to make sure it represents the bin size correctly
		bin_resolution = ((xi[1:]-xi[:-1]).mean() * positions.unit,(yi[1:]-yi[:-1]).mean() * positions.unit,(zi[1:]-zi[:-1]).mean() * positions.unit)

		#Perform smoothing if prompted
		if smooth is not None:

			#Fourier transform the density field
			fx,fy,fz = np.meshgrid(fftengine.fftfreq(density.shape[0]),fftengine.fftfreq(density.shape[1]),fftengine.rfftfreq(density.shape[2]),indexing="ij")
			density_ft = fftengine.rfftn(density)

			#Perform the smoothing
			density_ft *= np.exp(-0.5*((2.0*np.pi*smooth)**2)*(fx**2 + fy**2 + fz**2))

			#Go back in real space
			density = fftengine.irfftn(density_ft)


		#Return the density histogram, along with the bin resolution along each axis
		if save:
			self.density,self.resolution = density,bin_resolution

		return density,bin_resolution

	###################################################################################################################################################

	def cutPlaneGaussianGrid(self,normal=2,thickness=0.5*Mpc,center=7.0*Mpc,plane_resolution=4096,left_corner=None,thickness_resolution=1,smooth=1,kind="density",**kwargs):

		"""
		Cuts a density (or lensing potential) plane out of the snapshot by computing the particle number density on a slab and performing Gaussian smoothing; the plane coordinates are cartesian comoving

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

		:param kwargs: accepted keyword are: 'density_placeholder', a pre-allocated numpy array, with a RMA window opened on it; this facilitates the communication with different processors by using a single RMA window during the execution. 'l_squared' a pre-computed meshgrid of squared multipoles used for smoothing
		:type kwargs: dict.

		:returns: tuple(numpy 2D array with the density (or lensing potential),bin resolution along the axes, number of particles on the plane)

		"""

		#Sanity checks
		assert normal in range(3),"There are only 3 dimensions!"
		assert kind in ["density","potential"],"Specify density or potential plane!"
		assert type(thickness)==quantity.Quantity and thickness.unit.physical_type=="length"
		assert type(center)==quantity.Quantity and center.unit.physical_type=="length"

		#Redshift must be bigger than 0 or we cannot proceed
		if ("redshift" in self.header) and (self.header["redshift"]<=0.0):
			raise ValueError("The snapshot redshift must be >0 for the lensing density to be defined!")

		#Cosmological normalization factor
		cosmo_normalization = 1.5 * self.header["H0"]**2 * self.header["Om0"] / c**2

		#Direction of the plane
		plane_directions = [ d for d in range(3) if d!=normal ]

		#Get the particle positions if not available get
		if hasattr(self,"positions"):
			positions = self.positions
		else:
			positions = self.getPositions(first=self._first,last=self._last,save=False)

		assert hasattr(self,"weights")
		assert hasattr(self,"virial_radius")
		assert hasattr(self,"concentration")

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

		#Weights
		if self.weights is not None:
			weights = self.weights.astype(np.float32)
		else:
			weights = None

		#Virial radius
		if self.virial_radius is not None:
			assert weights is not None,"Particles have virial radiuses, you should specify their weight!"
			weights  = (weights * self._header["num_particles_total"] / ((len(binning[0]) - 1) * (len(binning[1]) - 1) * (len(binning[2]) - 1))).astype(np.float32)
			rv = self.virial_radius.to(positions.unit).value
		else:
			rv = None

		#Recompute resolution to make sure it represents the bin size correctly
		bin_resolution = [ (binning[n][1:]-binning[n][:-1]).mean() * positions.unit for n in (0,1,2) ]

		############################################################################################################
		#################################Longitudinal normalization factor##########################################
		#If the comoving distance is not provided in the header, the position along the normal direction is assumed#
		############################################################################################################

		if "comoving_distance" in self.header:
			
			#Constant time snapshots
			density_normalization = bin_resolution[normal] * self.header["comoving_distance"] / self.header["scale_factor"]
		
		else:

			#Light cone projection: use the lens center as the common comoving distance
			zlens = z_at_value(self.cosmology.comoving_distance,center)
			density_normalization = bin_resolution[normal] * center * (1.+zlens)

		#Now use gridding to compute the density along the slab
		assert positions.value.dtype==np.float32

		#Log
		if self.pool is not None:
			logplanes.debug("Task {0} began gridding procedure".format(self.pool.rank))
		else:
			logplanes.debug("Began gridding procedure")

		##########
		#Gridding#
		##########

		density = ext._nbody.grid3d_nfw(positions.value,tuple(binning),weights,rv,self.concentration)

		###################################################################################################################################

		#Log
		if self.pool is not None:
			logplanes.debug("Task {0} done with gridding procedure".format(self.pool.rank))
		else:
			logplanes.debug("Done with gridding procedure")

		if (self.pool is None) or (self.pool.is_master()):
			logstderr.debug("Done with gridding procedure: peak memory usage {0:.3f} (task)".format(peakMemory()))

		#Accumulate the density from the other processors
		if "density_placeholder" in kwargs.keys():

			density_projected = kwargs["density_placeholder"]

			#Safety assert
			assert density_projected.shape==(density.shape[plane_directions[0]],density.shape[plane_directions[1]])

			density_projected[:] = density.sum(normal)
			NumPartTask = density_projected.sum()

			if self.pool is not None:

				self.pool.comm.Barrier()

				#Log
				logplanes.debug("Task {0} collected {1:.3e} particles".format(self.pool.rank,NumPartTask))

				#Compute how many particles in total shoud be collected (for checking)
				NumPartTotalExpected = np.zeros(1,dtype=np.float32)
				self.pool.comm.Reduce(np.array([NumPartTask]),NumPartTotalExpected)

				#Log
				if self.pool.is_master():
					logplanes.debug("{0[0]:.3e} particles should be collected from tasks 0-{1}".format(NumPartTotalExpected,self.pool.size))
					logplanes.debug("Communicating density between tasks...")

				self.pool.accumulate()
		
		else:

			#Project along the normal direction
			density_projected = density.sum(normal)
			NumPartTask = density_projected.sum()
			
			if self.pool is not None:

				#Log
				logplanes.debug("Task {0} collected {1:.3e} particles".format(self.pool.rank,NumPartTask))
				
				self.pool.openWindow(density_projected)
				self.pool.accumulate()
				self.pool.closeWindow()

		#Safety barrier sync
		if self.pool is not None:
			self.pool.comm.Barrier()

		#Compute the number of particles on the plane
		NumPartTotal = density_projected.sum()

		#Log
		if (self.pool is not None) and self.pool.is_master():
			logplanes.debug("Received particles from all tasks: collected {0:.3e} particles".format(NumPartTotal))
			logstderr.debug("Received particles from all tasks: peak memory usage {0:.3f} (task)".format(peakMemory()))

		#If this task is not the master, we can return now
		if (self.pool is not None) and not(self.pool.is_master()):
			return (None,)*3

		#Normalize the density to the density fluctuation
		density_projected /= self._header["num_particles_total"]
		density_projected *= (self._header["box_size"]**3 / (bin_resolution[0]*bin_resolution[1]*bin_resolution[2])).decompose().value

		#################################################################################################################################
		######################################Ready to solve poisson equation via FFTs###################################################
		#################################################################################################################################

		bin_resolution.pop(normal)

		#If smoothing is enabled or potential calculations are needed, we need to FFT the density field
		if (smooth is not None) or kind=="potential":

			#Compute the multipoles
			if "l_squared" in kwargs.keys():

				l_squared = kwargs["l_squared"]
			
			else:
				lx,ly = np.meshgrid(fftengine.fftfreq(density_projected.shape[0]),fftengine.rfftfreq(density_projected.shape[1]),indexing="ij")
				l_squared = lx**2 + ly**2
				
				#Avoid dividing by 0
				l_squared[0,0] = 1.0


			#FFT the density field
			if (self.pool is None) or (self.pool.is_master()):
				logplanes.debug("Proceeding in density FFT operations...")

			density_ft = fftengine.rfftn(density_projected)

			#Zero out the zeroth frequency
			density_ft[0,0] = 0.0

			if kind=="potential":

				#Find out the comoving distance
				if "comoving_distance" in self.header:
					chi = self.header["comoving_distance"]
				else:
					chi = center

				#Solve the poisson equation
				density_ft *= -2.0 * (bin_resolution[0] * bin_resolution[1] / chi**2).decompose().value / (l_squared * ((2.0*np.pi)**2))

			if smooth is not None:
				#Perform the smoothing
				density_ft *= np.exp(-0.5*((2.0*np.pi*smooth)**2)*l_squared)

			#Revert the FFT
			lensing_potential = fftengine.irfftn(density_ft)

			if (self.pool is None) or (self.pool.is_master()):
				logplanes.debug("Done with density FFT operations...")
				logstderr.debug("Done with density FFT operations: peak memory usage {0:.3f} (task)".format(peakMemory()))


		else:

			lensing_potential = density_projected

		#Multiply by the normalization factors
		lensing_potential = lensing_potential * cosmo_normalization * density_normalization
		lensing_potential = lensing_potential.decompose()
		assert lensing_potential.unit.physical_type=="dimensionless"

		#Add units to lensing potential
		if kind=="potential":
			lensing_potential *= rad**2
		else:
			lensing_potential = lensing_potential.value

		#Return
		return lensing_potential,bin_resolution,NumPartTotal


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
			positions = self.positions
		else:
			positions = self.getPositions(save=False)

		assert hasattr(self,"weights")
		assert hasattr(self,"virial_radius")
		assert hasattr(self,"concentration")

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

		#Weights
		if self.weights is not None:
			weights = self.weights.astype(np.float32)
		else:
			weights = None

		#Compute the adaptive smoothing
		density = (3.0/np.pi)*ext._nbody.adaptive(positions.value,weights,rp,self.concentration,binning,center.to(positions.unit).value,plane_directions[0],plane_directions[1],normal,projectAll)

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
			lx,ly = np.meshgrid(fftengine.fftfreq(density.shape[0]),fftengine.rfftfreq(density.shape[1]),indexing="ij")
			l_squared = lx**2 + ly**2

			#Avoid dividing by 0
			l_squared[0,0] = 1.0

			#FFT the density field
			density_ft = fftengine.rfftn(density)

			#Zero out the zeroth frequency
			density_ft[0,0] = 0.0

			#Solve the poisson equation
			density_ft *= -2.0 * (bin_resolution[0] * bin_resolution[1] / self.header["comoving_distance"]**2).decompose().value / (l_squared * ((2.0*np.pi)**2))

			#Revert the FFT and return
			density = fftengine.irfftn(density_ft)
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

		assert hasattr(self,"weights")
		assert hasattr(self,"virial_radius")
		assert hasattr(self,"concentration")

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

		if self.virial_radius is not None:
			rv = self.virial_radius.to(positions.unit).value
		else:
			rv = None

		density = ext._nbody.grid3d(positions,tuple(binning),self.weights,rv,self.concentration)

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

				fx,fy,fz = np.meshgrid(fftengine.fftfreq(density.shape[0]),fftengine.fftfreq(density.shape[1]),fftengine.rfftfreq(density.shape[2]),indexing="ij")
				density_ft = fftengine.rfftn(density)
				density_ft *= np.exp(-0.5*((2.0*np.pi*smooth)**2)*(fx**2 + fy**2 + fz**2))
				density_ft[0,0] = 0.0
				density = fftengine.irfftn(density_ft)

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
			lx,ly = np.meshgrid(fftengine.fftfreq(density.shape[0]),fftengine.rfftfreq(density.shape[1]),indexing="ij")
			l_squared = lx**2 + ly**2

			#Avoid dividing by 0
			l_squared[0,0] = 1.0

			#Fourier transform the density field
			density_ft = fftengine.rfftn(density)

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
				density = fftengine.irfftn(density_ft)
			elif space=="fourier":
				density = density_ft
			else:
				raise ValueError("space must be real or fourier!")

		else:

			density -= density.sum() / reduce(mul,density.shape)
			if space=="fourier":
				density = fftengine.rfftn(density)

		#Return
		return (density*cosmo_normalization*density_normalization).decompose().value,bin_resolution,NumPartTotal


	#############################################################################################################################################

	def lensMaxSize(self):

		"""
		Computes the maximum observed size of a lens plane cut out of the current snapshot
	
		"""

		return ((self._header["box_size"] / self.cosmology.comoving_distance(self._header["redshift"])) * rad).to(deg)


	#############################################################################################################################################


	def powerSpectrum(self,k_edges,resolution=None,return_num_modes=False,density_placeholder=None):

		"""
		Computes the power spectrum of the relative density fluctuations in the snapshot at the wavenumbers specified by k_edges; a discrete particle number density is computed before hand to prepare the FFT grid

		:param k_edges: wavenumbers at which to compute the density power spectrum (must have units)
		:type k_edges: array.

		:param resolution: optional, fix the grid resolution to some value; to be passed to the massDensity method. If none this is computed automatically from the k_edges
		:type resolution: float with units, int. or None

		:param return_num_modes: if True returns the mode counting for each k bin as the last element in the return tuple
		:type return_num_modes: bool.

		:param density placeholder: if not None, it is used as a fixed memory chunk for MPI communications in the density calculations
		:type density_placeholder: array

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
			density,bin_resolution = self.massDensity(resolution=resolution,density_placeholder=density_placeholder) 
		else:
			assert resolution is None,"The spatial resolution is already specified in the attributes of this instance! Call massDensity() to modify!"
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
			logstderr.warning("Your grid resolution is too low to compute accurately the power on {0} (maximum recommended {1}, distortions might start to appear already at {2}): results might be inaccurate".format(k_edges.max(),k_max,k_max_recommended))

		#Perform the FFT
		density_ft = fftengine.rfftn(density)

		#Compute the azimuthal averages
		hits,power_spectrum = ext._topology.rfft3_azimuthal(density_ft,density_ft,kpixX.value,kpixY.value,kpixZ.value,k_edges.value)

		#Return the result (normalize the power so it corresponds to the one of the density fluctuations)
		k = 0.5*(k_edges[1:]+k_edges[:-1])
		return_tuple = (k,(power_spectrum/hits) * (bin_resolution[0] * bin_resolution[1] * bin_resolution[2])**2 / (self._header["box_size"]**3))

		if return_num_modes:
			return_tuple += (hits,)

		return return_tuple



	def __add__(self,rhs):

		"""
		Add two snapshots together: useful when the particle content is split between different files; all the positions and particle velocities are vstacked together

		"""

		merged_snapshot = self.__class__(None)
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



