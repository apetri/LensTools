"""

.. module:: flexion 
	:platform: Unix
	:synopsis: This module implements a set of operations which are usually performed on weak lensing flexion maps


.. moduleauthor:: Andrea Petri <apetri@phys.columbia.edu>
                ... and edited by Brij Patel <brp53@drexel.edu>

"""

from __future__ import division

from ..extern import _topology
from .convergence import ConvergenceMap

import numpy as np

#FFT engine
from ..utils.fft import NUMPYFFTPack
fftengine = NUMPYFFTPack()

#Units
from astropy.units import rad,arcsec,quantity

#I/O
from .io import loadFITS,saveFITS

try:
	import matplotlib
	import matplotlib.pyplot as plt
	matplotlib = matplotlib
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

		self._extra_attributes = kwargs.keys()
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


	#Multipole values in real FFT space
	def getEll(self):

		"""
		Get the values of the multipoles in real FFT space

		:returns: ell array with real FFT shape
		:rtype: array.

		"""

		ellx = fftengine.fftfreq(self.data.shape[1])*2.0*np.pi / self.resolution.to(u.rad).value
		elly = fftengine.rfftfreq(self.data.shape[1])*2.0*np.pi / self.resolution.to(u.rad).value
		return np.sqrt(ellx[:,None]**2 + elly[None,:]**2)

	###############################################################################################
	###############################################################################################

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


	
	def visualize(self,fig=None,ax=None,component_labels=("F1","F2"),colorbar=False,cmap="viridis",cbar_label=None,**kwargs):

		"""
		Visualize the flexion map; the kwargs are passed to imshow 

		"""

		if not matplotlib:
			raise ImportError("matplotlib is not installed, cannot visualize!")

		#Instantiate figure
		if (fig is None) or (ax is None):
			
			self.fig,self.ax = plt.subplots(1,self.data.shape[0],figsize=(16,8))

		else:

			self.fig = fig
			self.ax = ax

		#Build the color map
		if isinstance(cmap,matplotlib.colors.Colormap):
			cmap = cmap
		else:
			cmap = plt.get_cmap(cmap)

		#Plot the map
		if colorbar:

			for i in range(self.data.shape[0]):
				plt.colorbar(self.ax[i].imshow(self.data[i],origin="lower",interpolation="nearest",extent=[0,self.side_angle.value,0,self.side_angle.value],cmap=cmap,**kwargs),ax=self.ax[i])
				self.ax[i].grid(b=False)
		
		else:

			for i in range(self.data.shape[0]):
				self.ax[i].imshow(self.data[i],origin="lower",interpolation="nearest",extent=[0,self.side_angle.value,0,self.side_angle.value],cmap=cmap,**kwargs)
				self.ax[i].grid(b=False)

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


#############################################
##########FlexionMap class###################
#############################################

class FlexionMap(Spin1):

	"""
	A class that handles 2D flexion maps and allows to perform a set of operations on them

	"""

	#Construct flexion from convergence via KS ideology
	@classmethod
	def fromConvergence(cls,conv):

		"""
		Construct a flexion map from a ConvergenceMap instance using the Kaiser Squires ideology

		:param conv: input convergence map 
		:type conv: ConvergenceMap

		:returns: reconstructed flexion map
		:rtype: FlexionMap

		"""

		#Type check
		assert isinstance(conv,ConvergenceMap)

		#Multipoles
		lx = fftengine.rfftfreq(conv.data.shape[0])[None]
		ly = fftengine.fftfreq(conv.data.shape[0])[:,None]

		#FFT forward, rotation, FFT backwards
		conv_fft = fftengine.rfft2(conv.data)
		F1 = fftengine.irfft2(1j*lx*conv_fft)
		F2 = fftengine.irfft2(1j*ly*conv_fft)

		#Return
		kwargs = dict((k,getattr(conv,k)) for k in conv._extra_attributes)
		return cls(np.array([F1,F2]),conv.side_angle,**kwargs)

	#Construct convergence map with KS
	def convergence(self):
		
		"""
		Reconstructs the convergence from flexion using Kaiser Squires ideology

		:returns: new ConvergenceMap instance 

		"""

     #Perform Fourier transforms
		ft_F1 = fftengine.rfft2(self.data[0])
		ft_F2 = fftengine.rfft2(self.data[1])

		#Compute frequencies
		lx = fftengine.rfftfreq(ft_F1.shape[0])
		ly = fftengine.fftfreq(ft_F1.shape[0])

		#Safety check
		assert len(lx)==ft_F1.shape[1]
		assert len(ly)==ft_F1.shape[0]

		l_squared = lx[np.newaxis,:]**2 + ly[:,np.newaxis]**2
		l_squared[0,0] = 1.0

		#Compute Fourier Transform of the convergence
		ft_conv = -1j*(ft_F1*lx[np.newaxis,:] + ft_F2*ly[:,np.newaxis])/l_squared
		ft_conv[0,0] = 0.0

		assert ft_conv.shape == ft_F1.shape

		#Invert the Fourier transform to go back to real space to get real convergence
		conv = fftengine.irfft2(ft_conv)

		#Return the ConvergenceMap instance
		kwargs = dict((k,getattr(self,k)) for k in self._extra_attributes)
		return ConvergenceMap(conv,self.side_angle,**kwargs)