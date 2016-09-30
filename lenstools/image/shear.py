"""

.. module:: shear 
	:platform: Unix
	:synopsis: This module implements a set of operations which are usually performed on weak lensing shear maps


.. moduleauthor:: Andrea Petri <apetri@phys.columbia.edu>


"""

from __future__ import division

from ..extern import _topology
from .convergence import ConvergenceMap

import numpy as np

#FFT engine
from ..utils.fft import NUMPYFFTPack
fftengine = NUMPYFFTPack()

#Units
from astropy.units import deg,rad,arcsec,quantity

#I/O
from .io import loadFITS,saveFITS

try:
	import matplotlib.pyplot as plt
	from matplotlib.colors import LogNorm
	matplotlib = True
except ImportError:
	matplotlib = False


##########################################
########Spin1 class#######################
##########################################

class Spin1(object):


	def __init__(self,data,angle,**kwargs):

		#Sanity check
		assert angle.unit.physical_type in ["angle","length"]
		assert data.shape[1]==data.shape[2],"The map must be a square!!"

		self.data = data
		self.side_angle = angle
		self.resolution = self.side_angle / self.data.shape[1]

		if self.side_angle.unit.physical_type=="angle":
			self.resolution = self.resolution.to(arcsec)
			self.lmin = 2.0*np.pi/self.side_angle.to(rad).value
			self.lmax = np.sqrt(2)*np.pi/self.resolution.to(rad).value

		self._extra_args = kwargs.keys()
		for key in kwargs:
			setattr(self,key,kwargs[key])

	@property
	def info(self):

		"""
		Displays some of the information stored in the map (mainly resolution)

		"""

		print("Pixels on a side: {0}".format(self.data.shape[1]))
		print("Pixel size: {0}".format(self.resolution))
		print("Total angular size: {0}".format(self.side_angle))
		print("lmin={0:.1e} ; lmax={1:.1e}".format(self.lmin,self.lmax))
	

	@classmethod
	def load(cls,filename,format=None,**kwargs):
		
		"""
		
		This class method allows to read the map from a data file, in various formats

		:param filename: name of the file in which the map is saved
		:type filename: str. 
		
		:param format: the format of the file in which the map is saved (can be a callable too); if None, it's detected automatically from the filename
		:type format: str. or callable

		:param kwargs: the keyword arguments are passed to the format (if callable)
		:type kwargs: dict.

		:returns: Spin1 instance with the loaded map
		
		"""

		if format is None:
			
			extension = filename.split(".")[-1]
			if extension in ["fit","fits"]:
				format="fits"
			else:
				raise IOError("File format not recognized from extension '{0}', please specify it manually".format(extension))


		if format=="fits":
			return loadFITS(cls,filename)

		else:

			angle,data = format(filename,**kwargs)
			return cls(data,angle)


	def save(self,filename,format=None,double_precision=False):

		"""
		Saves the map to an external file, of which the format can be specified (only fits implemented so far)

		:param filename: name of the file on which to save the plane
		:type filename: str.

		:param format: format of the file, only FITS implemented so far; if None, it's detected automatically from the filename
		:type format: str.

		:param double_precision: if True saves the Plane in double precision
		:type double_precision: bool.

		"""

		if format is None:
			
			extension = filename.split(".")[-1]
			if extension in ["fit","fits"]:
				format="fits"
			else:
				raise IOError("File format not recognized from extension '{0}', please specify it manually".format(extension))


		if format=="fits":
			saveFITS(self,filename,double_precision)

		else:
			raise ValueError("Format {0} not implemented yet!!".format(format))


	def setAngularUnits(self,unit):

		"""
		Convert the angular units of the map to the desired unit

		:param unit: astropy unit instance to which to perform the conversion
		:type unit: astropy units 
		
		"""

		#Sanity check
		assert unit.physical_type=="angle"
		self.side_angle = self.side_angle.to(unit)


	def gradient(self,x=None,y=None):

		"""
		Computes the gradient of the components of the spin1 field at each point

		:param x: optional, x positions at which to evaluate the gradient
		:type x: array with units

		:param y: optional, y positions at which to evaluate the gradient
		:type y: array with units

		:returns: the gradient of the spin1 field in array form, of shape (4,:,:) where the four components are, respectively, 1x,1y,2x,2y; the units for the finite difference are pixels

		"""

		if self.data.shape[0] > 2:
			raise ValueError("Gradients are nor defined yet for spin>1 fields!!")

		if (x is not None) and (y is not None):

			assert x.shape==y.shape,"x and y must have the same shape!"

			#x coordinates
			if type(x)==quantity.Quantity:
			
				assert x.unit.physical_type=="angle"
				j = np.mod(((x / self.resolution).decompose().value).astype(np.int32),self.data.shape[1])

			else:

				j = np.mod((x / self.resolution.to(rad).value).astype(np.int32),self.data.shape[1])	

			#y coordinates
			if type(y)==quantity.Quantity:
			
				assert y.unit.physical_type=="angle"
				i = np.mod(((y / self.resolution).decompose().value).astype(np.int32),self.data.shape[1])

			else:

				i = np.mod((y / self.resolution.to(rad).value).astype(np.int32),self.data.shape[1])

		else:
			i = None
			j = None

		#Call the C backend
		grad1x,grad1y = _topology.gradient(self.data[0],j,i)
		grad2x,grad2y = _topology.gradient(self.data[1],j,i)

		#Return
		if (x is not None) and (y is not None):
			return np.array([grad1x.reshape(x.shape),grad1y.reshape(y.shape),grad2x.reshape(x.shape),grad2y.reshape(y.shape)])
		else:
			return np.array([grad1x,grad1y,grad2x,grad2y])


	def getValues(self,x,y):

		"""
		Extract the map values at the requested (x,y) positions; this is implemented using the numpy fast indexing routines, so the formats of x and y must follow the numpy advanced indexing rules. Periodic boundary conditions are enforced

		:param x: x coordinates at which to extract the map values (if unitless these are interpreted as radians)
		:type x: numpy array or quantity 

		:param y: y coordinates at which to extract the map values (if unitless these are interpreted as radians)
		:type y: numpy array or quantity 

		:returns: numpy array with the map values at the specified positions, with shape (N,shape x) where N is the number of components of the map field

		:raises: IndexError if the formats of x and y are not the proper ones

		"""

		assert isinstance(x,np.ndarray) and isinstance(y,np.ndarray)

		#x coordinates
		if type(x)==quantity.Quantity:
			
			assert x.unit.physical_type=="angle"
			j = np.mod(((x / self.resolution).decompose().value).astype(np.int32),self.data.shape[2])

		else:

			j = np.mod((x / self.resolution.to(rad).value).astype(np.int32),self.data.shape[2])	

		#y coordinates
		if type(y)==quantity.Quantity:
			
			assert y.unit.physical_type=="angle"
			i = np.mod(((y / self.resolution).decompose().value).astype(np.int32),self.data.shape[1])

		else:

			i = np.mod((y / self.resolution.to(rad).value).astype(np.int32),self.data.shape[1])

		#Return the map values at the specified coordinates
		return self.data[:,i,j]


	
	def visualize(self,fig=None,ax=None,component_labels=(r"$\gamma_1$",r"$\gamma_2$"),colorbar=False,cmap="jet",cbar_label=None,**kwargs):

		"""
		Visualize the shear map; the kwargs are passed to imshow 

		"""

		if not matplotlib:
			raise ImportError("matplotlib is not installed, cannot visualize!")

		#Instantiate figure
		if (fig is None) or (ax is None):
			
			self.fig,self.ax = plt.subplots(1,self.data.shape[0],figsize=(16,8))

		else:

			self.fig = fig
			self.ax = ax

		#Plot the map
		if colorbar:

			for i in range(self.data.shape[0]):
				plt.colorbar(self.ax[i].imshow(self.data[i],origin="lower",interpolation="nearest",extent=[0,self.side_angle.value,0,self.side_angle.value],**kwargs),ax=self.ax[i])
		
		else:

			for i in range(self.data.shape[0]):
				self.ax[i].imshow(self.data[i],origin="lower",interpolation="nearest",extent=[0,self.side_angle.value,0,self.side_angle.value],cmap=plt.get_cmap(cmap),**kwargs)

		#Axes labels
		for i in range(self.data.shape[0]):

			self.ax[i].set_xlabel(r"$x$({0})".format(self.side_angle.unit.to_string()),fontsize=18)
			self.ax[i].set_ylabel(r"$y$({0})".format(self.side_angle.unit.to_string()),fontsize=18)
			self.ax[i].set_title(component_labels[i],fontsize=18)

	
	def savefig(self,filename):

		"""
		Saves the map visualization to an external file

		:param filename: name of the file on which to save the map
		:type filename: str.

		"""

		self.fig.savefig(filename)



##########################################
########Spin2 class#######################
##########################################

class Spin2(Spin1):

	@classmethod
	def fromEBmodes(cls,fourier_E,fourier_B,angle=3.14*deg):

		"""
		This class method allows to build a shear map specifying its E and B mode components

		:param fourier_E: E mode of the shear map in fourier space
		:type fourier_E: numpy 2D array, must be of type np.complex128 and must have a shape that is appropriate for a real fourier transform, i.e. (N,N/2 + 1); N should be a power of 2

		:param fourier_B: B mode of the shear map in fourier space
		:type fourier_B: numpy 2D array, must be of type np.complex128 and must have a shape that is appropriate for a real fourier transform, i.e. (N,N/2 + 1); N should be a power of 2

		:param angle: Side angle of the real space map in degrees
		:type angle: float.

		:returns: the corresponding ShearMap instance

		:raises: AssertionErrors for inappropriate inputs

		"""

		assert fourier_E.dtype == np.complex128 and fourier_B.dtype == np.complex128
		assert fourier_E.shape[1] == fourier_E.shape[0]/2 + 1
		assert fourier_B.shape[1] == fourier_B.shape[0]/2 + 1
		assert fourier_E.shape == fourier_B.shape

		#Compute frequencies
		lx = fftengine.rfftfreq(fourier_E.shape[0])
		ly = fftengine.fftfreq(fourier_E.shape[0])

		#Safety check
		assert len(lx)==fourier_E.shape[1]
		assert len(ly)==fourier_E.shape[0]

		#Compute sines and cosines of rotation angles
		l_squared = lx[np.newaxis,:]**2 + ly[:,np.newaxis]**2
		l_squared[0,0] = 1.0

		sin_2_phi = 2.0 * lx[np.newaxis,:] * ly[:,np.newaxis] / l_squared
		cos_2_phi = (lx[np.newaxis,:]**2 - ly[:,np.newaxis]**2) / l_squared

		sin_2_phi[0,0] = 0.0
		cos_2_phi[0,0] = 0.0

		#Invert E/B modes and find the components of the shear
		ft_data1 = cos_2_phi * fourier_E - sin_2_phi * fourier_B
		ft_data2 = sin_2_phi * fourier_E + cos_2_phi * fourier_B

		#Invert Fourier transforms
		data1 = fftengine.irfft2(ft_data1)
		data2 = fftengine.irfft2(ft_data2)

		#Instantiate new shear map class
		new = cls(np.array([data1,data2]),angle)
		setattr(new,"fourier_E",fourier_E)
		setattr(new,"fourier_B",fourier_B)

		return new


	def sticks(self,fig=None,ax=None,pixel_step=10,multiplier=1.0):

		"""
		Draw the ellipticity map using the shear components

		:param ax: ax on which to draw the ellipticity field
		:type ax: matplotlib ax object

		:param pixel_step: One arrow will be drawn every pixel_step pixels to avoid arrow overplotting
		:type pixel_step: int.

		:param multiplier: Multiplies the stick length by a factor
		:type multiplier: float.

		:returns: ax -- the matplotlib ax object on which the stick field was drawn

		>>> import matplotlib.pyplot as plt
		>>> test = ShearMap.load("shear.fit",loader=load_fits_default_shear)
		>>> fig,ax = plt.subplots()
		>>> test.sticks(ax,pixel_step=50)

		"""

		if not(matplotlib):
			raise ImportError("matplotlib is not installed, cannot visualize!")

		#Instantiate fig,ax objects
		if (fig is None) or (ax is None):
			self.fig,self.ax = plt.subplots()
		else:
			self.fig = fig
			self.ax = ax

		x,y = np.meshgrid(np.arange(0,self.data.shape[2],pixel_step),np.arange(0,self.data.shape[1],pixel_step))

		#Translate shear components into sines and cosines
		cos_2_phi = self.data[0] / np.sqrt(self.data[0]**2 + self.data[1]**2)
		sin_2_phi = self.data[1] / np.sqrt(self.data[0]**2 + self.data[1]**2)

		#Compute stick directions
		cos_phi = np.sqrt(0.5*(1.0 + cos_2_phi)) * np.sign(sin_2_phi)
		sin_phi = np.sqrt(0.5*(1.0 - cos_2_phi))

		#Fix ambiguity when sin_2_phi = 0
		cos_phi[sin_2_phi==0] = np.sqrt(0.5*(1.0 + cos_2_phi[sin_2_phi==0]))

		#Draw map using matplotlib quiver
		self.ax.quiver(x*self.side_angle.value/self.data.shape[2],y*self.side_angle.value/self.data.shape[1],cos_phi[x,y],sin_phi[x,y],headwidth=0,units="height",scale=x.shape[0]/multiplier)

		#Axes labels
		self.ax.set_xlabel(r"$x$({0})".format(self.side_angle.unit.to_string()))
		self.ax.set_ylabel(r"$y$({0})".format(self.side_angle.unit.to_string()))



	def fourierEB(self):

		"""
		Computes E and B modes of the shear map in Fourier space

		:returns: (E,B) map
		:rtype: tuple.

		"""

		#Perform Fourier transforms
		ft_data1 = fftengine.rfft2(self.data[0])
		ft_data2 = fftengine.rfft2(self.data[1])

		#Compute frequencies
		lx = fftengine.rfftfreq(ft_data1.shape[0])
		ly = fftengine.fftfreq(ft_data1.shape[0])

		#Safety check
		assert len(lx)==ft_data1.shape[1]
		assert len(ly)==ft_data1.shape[0]

		#Compute sines and cosines of rotation angles
		l_squared = lx[np.newaxis,:]**2 + ly[:,np.newaxis]**2
		l_squared[0,0] = 1.0

		sin_2_phi = 2.0 * lx[np.newaxis,:] * ly[:,np.newaxis] / l_squared
		cos_2_phi = (lx[np.newaxis,:]**2 - ly[:,np.newaxis]**2) / l_squared

		#Compute E and B components
		ft_E = cos_2_phi * ft_data1 + sin_2_phi * ft_data2
		ft_B = -1.0 * sin_2_phi * ft_data1 + cos_2_phi * ft_data2

		ft_E[0,0] = 0.0
		ft_B[0,0] = 0.0

		assert ft_E.shape == ft_B.shape
		assert ft_E.shape == ft_data1.shape

		#Return 
		return ft_E,ft_B


	def eb_power_spectrum(self,l_edges):

		"""
		Decomposes the shear map into its E and B modes components and returns the respective power spectral densities at the specified multipole moments

		:param l_edges: Multipole bin edges
		:type l_edges: array

		:returns: (l -- array,P_EE,P_BB,P_EB -- arrays) = (multipole moments, EE,BB power spectra and EB cross power)
		:rtype: tuple.

		>>> test_map = ShearMap.load("shear.fit",format=load_fits_default_shear)
		>>> l_edges = np.arange(300.0,5000.0,200.0)
		>>> l,EE,BB,EB = test_map.eb_power_spectrum(l_edges)

		"""

		#Compute the E,B modes in Fourier space
		ft_E,ft_B = self.fourierEB()

		#Compute and return power spectra
		l = 0.5*(l_edges[:-1] + l_edges[1:])
		P_ee = _topology.rfft2_azimuthal(ft_E,ft_E,self.side_angle.to(deg).value,l_edges)
		P_bb = _topology.rfft2_azimuthal(ft_B,ft_B,self.side_angle.to(deg).value,l_edges)
		P_eb = _topology.rfft2_azimuthal(ft_E,ft_B,self.side_angle.to(deg).value,l_edges)

		#Return to user
		return l,P_ee,P_bb,P_eb

	def decompose(self,l_edges):
		return self.eb_power_spectrum(l_edges)

	def visualizeComponents(self,fig,ax,components="EE,BB,EB",region=(200,9000,-9000,9000)):

		"""
		Plots the full 2D E and B mode power spectrum (useful to test statistical isotropicity)

		:param fig: figure on which to draw the ellipticity field
		:type fig: matplotlibfigure object

		:param ax: ax on which to draw the ellipticity field
		:type ax: matplotlib ax object or array of ax objects, can be None in which case new ax are created

		:param components: string that contains the components to plot; the format is a sequence of {EE,BB,EB} separated by commas 
		:type components:

		:param region: selects the multipole region to visualize
		:type region: tuple (lx_min,lx_max,ly_min,ly_max)

		"""

		if not(matplotlib):
			raise ImportError("matplotlib is not installed, cannot visualize!")

		#Set value for frequency pixelization
		lpix = 360.0/self.side_angle.to(deg).value

		#First parse the components to plot from the components string
		component_list = components.split(",")
		if not len(component_list):
			return None

		for component in component_list:
			assert component=="EE" or component=="BB" or component=="EB", "Each of the components should be one of {EE,BB,EB}"

		#Compute EB modes in Fourier space
		ft_E,ft_B = self.fourierEB()

		#Instantiate figure
		if (fig is None) or (ax is None):
			
			self.fig,self.ax = plt.subplots(1,len(component_list))

		else:

			assert len(component_list)==np.size(ax), "You should specify a plotting ax for each component!"
			self.fig = fig
			self.ax = ax


		#Do the plotting
		for n,component in enumerate(component_list):

			if len(component_list)==1:
				plot_ax = self.ax
			else:
				plot_ax = self.ax[n]

			#Select the numpy array of the appropriate component
			if component=="EE":
				mode = np.abs(ft_E)**2 
			elif component=="BB":
				mode = np.abs(ft_B)**2
			elif component=="EB":
				mode = np.abs((ft_E * ft_B.conjugate()).real)

			#Roll it to get the frequencies in the right order
			mode = np.roll(mode,mode.shape[0]//2,axis=0)  * (self.side_angle.to(rad).value)**2/(ft_E.shape[0]**4)

			#Plot the components with the right frequencies on the axes
			plot_cbar = plot_ax.imshow(mode,origin="lower",interpolation="nearest",norm=LogNorm(),extent=[0,lpix*mode.shape[1],-lpix*mode.shape[0]/2,lpix*mode.shape[0]/2])
			plot_ax.set_xlim(region[0],region[1])
			plot_ax.set_ylim(region[2],region[3])

			#Set labels
			plot_ax.set_xlabel(r"$l_x$")
			plot_ax.set_ylabel(r"$l_y$")
			plot_ax.set_title(r"${0}$".format(component))

			#Set colorbar
			plt.colorbar(plot_cbar,ax=plot_ax)

			#Set tick size
			plot_ax.tick_params(labelsize='small')
			plot_ax.xaxis.get_major_formatter().set_powerlimits((0, 1))
			plot_ax.yaxis.get_major_formatter().set_powerlimits((0, 1))
			plot_ax.ticklabel_format(style='sci')


#############################################
###########ShearMap class####################
#############################################

class ShearMap(Spin2):

	"""
	A class that handles 2D shear maps and allows to perform a set of operations on them

	>>> from lenstools import ShearMap
	
	>>> test = ShearMap.load("shear.fit",format=lenstools.defaults.load_fits_default_shear)
	>>> test.side_angle
	1.95
	>>> test.data
	#The actual map values

	"""

	def convergence(self):
		
		"""
		Reconstructs the convergence from the E component of the shear

		:returns: new ConvergenceMap instance 

		"""

		#Compute EB modes in fourier space
		ft_E,ft_B = self.fourierEB()

		#Invert the Fourier transform to go back to real space
		conv = fftengine.irfft2(ft_E)

		#Return the ConvergenceMap instance
		return ConvergenceMap(conv,self.side_angle)


