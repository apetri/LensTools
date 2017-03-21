"""

.. module:: convergence
	:platform: Unix
	:synopsis: This module implements the tools to compute topological statistics on 2D convergence maps: measure the power spectrum, counting peaks, measure the minkowski functionals; handling of masks is also provided


.. moduleauthor:: Andrea Petri <apetri@phys.columbia.edu>


"""

from __future__ import division

from operator import mul
from functools import reduce
import numbers

from ..extern import _topology

import numpy as np
import scipy.special as sp

#FFT engine
from ..utils.fft import NUMPYFFTPack
fftengine = NUMPYFFTPack()

#Hankel transform
from ..utils import fht

from scipy.ndimage import filters

#Units
import astropy.units as u

#I/O
from .io import loadFITS,saveFITS,loadNPZ,saveNPZ

#CMB lensing
from .cmblens import QuickLens as Lens

#Plotting
try:
	import matplotlib
	import matplotlib.pyplot as plt
	from mpl_toolkits.axes_grid1 import make_axes_locatable
	matplotlib = matplotlib
except ImportError:
	matplotlib = False


################################################
########Spin0 class#############################
################################################

class Spin0(object):

	def __init__(self,data,angle,masked=False,**kwargs):

		#Sanity check
		assert angle.unit.physical_type in ["angle","length"]

		if not(data.dtype==np.complex):
			assert data.shape[0]==data.shape[1],"The map must be a square!!"

		#Convert to double precision for calculation accuracy
		if data.dtype in (np.float,np.complex):
			self.data = data
		else:
			self.data = data.astype(np.float)
			
		self.side_angle = angle
		self.resolution = self.side_angle / self.data.shape[0]

		if self.side_angle.unit.physical_type=="angle":
			
			self.resolution = self.resolution.to(u.arcsec)
			self.lmin = 2.0*np.pi/self.side_angle.to(u.rad).value
			self.lmax = np.sqrt(2)*np.pi/self.resolution.to(u.rad).value
		
		self._masked = masked

		#Extra keyword arguments
		self._extra_attributes = kwargs.keys()
		for key in kwargs:
			setattr(self,key,kwargs[key])

	@property
	def info(self):

		"""
		Displays some of the information stored in the map (mainly resolution)

		"""

		print("Pixels on a side: {0}".format(self.data.shape[0]))
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

		ellx = fftengine.fftfreq(self.data.shape[0])*2.0*np.pi / self.resolution.to(u.rad).value
		elly = fftengine.rfftfreq(self.data.shape[0])*2.0*np.pi / self.resolution.to(u.rad).value
		return np.sqrt(ellx[:,None]**2 + elly[None,:]**2)


	####################################################################################################
	####################################################################################################

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

		:returns: Spin0 instance with the loaded map

		"""

		if format is None:
			
			extension = filename.split(".")[-1]
			if extension in ["fit","fits"]:
				format="fits"
			elif extension in ["npy","npz"]:
				format="npz"
			else:
				raise IOError("File format not recognized from extension '{0}', please specify it manually".format(extension))

		if format=="fits":
			return loadFITS(cls,filename)
		elif format=="npz":
			return loadNPZ(cls,filename)
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
			elif extension in ["npy","npz"]:
				format="npz"
			else:
				raise IOError("File format not recognized from extension '{0}', please specify it manually".format(extension))

		
		if format=="fits":
			saveFITS(self,filename,double_precision)
		elif format=="npz":
			saveNPZ(self,filename)
		else:
			raise ValueError("Format {0} not implemented yet!!".format(format))

	####################################################################################################
	####################################################################################################

	def setAngularUnits(self,unit):

		"""
		Convert the angular units of the map to the desired unit

		:param unit: astropy unit instance to which to perform the conversion
		:type unit: astropy units 
		
		"""

		#Sanity check
		assert unit.physical_type=="angle"
		self.side_angle = self.side_angle.to(unit)

	####################################################################################################
	####################################################################################################

	def mean(self):

		"""
		Computes the mean value of the pixels, taking into account eventual masking

		:returns: float.

		"""

		if not self._masked:
			
			return self.data.mean()
		
		else:
			
			if not hasattr(self,"_full_mask"):
				self.maskBoundaries()
			
			return self.data[self._full_mask].mean()


	def std(self):

		"""
		Computes the standard deviation of the pixels, taking into account eventual masking

		:returns: float.

		"""

		if not self._masked:

			return self.data.std()

		else:

			if not hasattr(self,"_full_mask"):
				self.maskBoundaries()

			return self.data[self._full_mask].std()


	def getValues(self,x,y):

		"""
		Extract the map values at the requested (x,y) positions; this is implemented using the numpy fast indexing routines, so the formats of x and y must follow the numpy advanced indexing rules. Periodic boundary conditions are enforced

		:param x: x coordinates at which to extract the map values (if unitless these are interpreted as radians)
		:type x: numpy array or quantity 

		:param y: y coordinates at which to extract the map values (if unitless these are interpreted as radians)
		:type y: numpy array or quantity 

		:returns: numpy array with the map values at the specified positions, with the same shape as x and y

		:raises: IndexError if the formats of x and y are not the proper ones

		"""

		assert isinstance(x,np.ndarray) and isinstance(y,np.ndarray)

		#x coordinates
		if type(x)==u.quantity.Quantity:
			
			assert x.unit.physical_type==self.side_angle.unit.physical_type
			j = np.mod(((x / self.resolution).decompose().value).astype(np.int32),self.data.shape[1])

		else:

			j = np.mod((x / self.resolution.to(u.rad).value).astype(np.int32),self.data.shape[1])	

		#y coordinates
		if type(y)==u.quantity.Quantity:
			
			assert y.unit.physical_type==self.side_angle.unit.physical_type
			i = np.mod(((y / self.resolution).decompose().value).astype(np.int32),self.data.shape[0])

		else:

			i = np.mod((y / self.resolution.to(u.rad).value).astype(np.int32),self.data.shape[0])

		#Return the map values at the specified coordinates
		return self.data[i,j]


	def cutRegion(self,extent):

		"""
		Cut a subset of the map, with the same resolution (warning! tested on square cuts only!)

		:param extent: region in the map to select (xmin,xmax,ymin,ymax), must have units
		:type extent: array with units

		:returns: new ConvergenceMap instance that encloses the selected region

		"""

		xmin,xmax,ymin,ymax = extent
		assert (xmax-xmin)==(ymax-ymin),"Only square cuts implemented so far!"

		x = np.arange(xmin.value,xmax.value,self.resolution.to(extent.unit).value)
		y = np.arange(ymin.value,ymax.value,self.resolution.to(extent.unit).value)

		#Initialize the meshgrid
		xx,yy = np.meshgrid(x,y) * extent.unit

		return self.__class__(data=self.getValues(xx,yy),angle=(self.resolution*xx.shape[0]).to(extent.unit),_extent=extent)



	def visualize(self,fig=None,ax=None,colorbar=False,cmap="viridis",cbar_label=None,**kwargs):

		"""
		Visualize the convergence map; the kwargs are passed to imshow 

		"""

		if not matplotlib:
			raise ImportError("matplotlib is not installed, cannot visualize!")

		#Instantiate figure
		if (fig is None) or (ax is None):
			
			self.fig,self.ax = plt.subplots()

		else:

			self.fig = fig
			self.ax = ax

		#Plot the map
		if hasattr(self,"_extent"):
			extent = self._extent.to(self.side_angle.unit).value
		else:
			extent = [0,self.side_angle.value,0,self.side_angle.value]

		#Build the color map
		if isinstance(cmap,matplotlib.colors.Colormap):
			cmap = cmap
		else:
			cmap = plt.get_cmap(cmap)

		#Visualize
		ax0 = self.ax.imshow(self.data,origin="lower",interpolation="nearest",extent=extent,cmap=cmap,**kwargs)
		self.ax.grid(b=False)
		self.ax.set_xlabel(r"$x$({0})".format(self.side_angle.unit.to_string()),fontsize=18)
		self.ax.set_ylabel(r"$y$({0})".format(self.side_angle.unit.to_string()),fontsize=18)

		#Add a colorbar maybe
		if colorbar:

			divider = make_axes_locatable(self.ax)
			cax = divider.append_axes("right",size="5%",pad=0.05)
			cbar = plt.colorbar(ax0,cax=cax)

			if cbar_label is not None:
				cbar.set_label(r"{0}".format(cbar_label),size=18)

	def savefig(self,filename):

		"""
		Saves the map visualization to an external file

		:param filename: name of the file on which to save the map
		:type filename: str.

		"""

		self.fig.savefig(filename)

	def mask(self,mask_profile,inplace=False):

		"""
		Applies a mask to the convergence map: all masked pixels are given a nan value because they cannot be used in the calculations

		:param mask_profile: profile of the mask, must be an array of 1 byte intergers that are either 0 (if the pixel is masked) or 1 (if the pixel is not masked). Must be of the same shape as the original map
		:type mask_profile: array. or ConvergenceMap instance

		:param inplace: if True the masking is performed in place and the original map is lost, otherwise a new instance of ConvergenceMap with the masked map is returned

		:returns: the masked convergence map if inplace is False, otherwise a float corresponding to the masked fraction of the map

		"""

		if inplace:
			new_map = self
		else:
			new_map = self.__class__(self.data.copy(),self.side_angle) 

		if isinstance(mask_profile,np.ndarray):

			assert mask_profile.dtype == np.int8
			assert len(mask_profile[mask_profile==0]) + len(mask_profile[mask_profile==1]) == reduce(mul,mask_profile.shape),"The mask must be made of 0 and 1 only!"
			assert mask_profile.shape == self.data.shape

			new_map.data[mask_profile==0] = np.nan

		elif isinstance(mask_profile,self.__class__):

			assert mask_profile.side_angle == self.side_angle
			assert len(mask_profile.data[mask_profile.data==0]) + len(mask_profile.data[mask_profile.data==1]) == reduce(mul,mask_profile.data.shape),"The mask must be made of 0 and 1 only!"
			assert mask_profile.data.shape == self.data.shape

			new_map.data[mask_profile.data==0] = np.nan

		else: 

			raise TypeError("Mask type not supported")

		new_map._masked = True
		new_map._mask = ~np.isnan(new_map.data)
		new_map._masked_fraction = 1.0 - new_map._mask.sum() / reduce(mul,new_map.data.shape)

		#Recompute gradients
		if (hasattr(new_map,"gradient_x") or hasattr(new_map,"gradient_y")):
			new_map.gradient()

		if (hasattr(new_map,"hessian_xx") or hasattr(new_map,"hessian_yy") or hasattr(new_map,"hessian_xy")):
			new_map.hessian()

		#Return
		if inplace:
			return new_map._masked_fraction
		else:
			return new_map 

	def maskBoundaries(self):

		"""
		Computes the mask boundaries defined in the following way: a boundary is a region where the convergence value is defined, but the gradients are not defined.

		:returns: float. : perimeter/area ratio of the mask boundaries
		
		"""

		if not self._masked:
			print("The map is not masked!!")
			return None

		#First check that the instance has the gradient and hessian attributes; if not, compute them
		if not (hasattr(self,"gradient_x") and hasattr(self,"gradient_y")):
			self.gradient()

		if not (hasattr(self,"hessian_xx") and hasattr(self,"hessian_yy") and hasattr(self,"hessian_xy")):
			self.hessian()

		#Check where gradient starts to have problems
		nan_gradient_pixels = np.isnan(self.gradient_x) + np.isnan(self.gradient_y)
		self._gradient_boundary = ~self._mask ^ nan_gradient_pixels

		#Check where the hessian has alsp problems
		nan_gradient_pixels = nan_gradient_pixels + np.isnan(self.hessian_xx) + np.isnan(self.hessian_yy) + np.isnan(self.hessian_xy)
		self._hessian_boundary = ~self._mask ^ nan_gradient_pixels

		#Create attribute that holds the full mask (including gradients)
		self._full_mask = self._mask * (~nan_gradient_pixels)
		assert self._full_mask.sum() < self._mask.sum()

		#Compute perimeter/area of the mask
		perimeter_area = self._hessian_boundary.sum() / (~self._full_mask).sum()

		#Return
		return perimeter_area

	@property
	def maskedFraction(self):

		"""
		Computes the masked fraction of the map

		:returns: float, the masked fraction of the map

		"""

		if not self._masked:
			return 0.0
		else:
			return self._masked_fraction

	@property
	def boundary(self):

		"""
		Computes the boundaries of the masked regions, defined as the regions in which the convergence is still well defined but the first and second derivatives are not

		:returns: array of bool. of the same shape as the map, with True values along the boundaries
		
		"""

		if not hasattr(self,"_hessian_boundary"):
			self.maskBoundaries()

		return self._hessian_boundary
			

	def gradient(self,x=None,y=None,save=True):
		
		"""
		Computes the gradient of the map and sets the gradient_x,gradient_y attributes accordingly

		:param x: optional, x positions at which to evaluate the gradient
		:type x: array with units

		:param y: optional, y positions at which to evaluate the gradient
		:type y: array with units

		:param save: if True saves the gradient as attrubutes
		:type save: bool.

		:returns: tuple -- (gradient_x,gradient_y)

		>>> test_map = ConvergenceMap.load("map.fit")
		>>> gx,gy = test_map.gradient()

		"""

		if (x is not None) and (y is not None):

			assert x.shape==y.shape,"x and y must have the same shape!"

			#x coordinates
			if type(x)==u.quantity.Quantity:
			
				assert x.unit.physical_type==self.side_angle.unit.physical_type
				j = np.mod(((x / self.resolution).decompose().value).astype(np.int32),self.data.shape[1])

			else:

				j = np.mod((x / self.resolution.to(u.rad).value).astype(np.int32),self.data.shape[1])	

			#y coordinates
			if type(y)==u.quantity.Quantity:
			
				assert y.unit.physical_type==self.side_angle.unit.physical_type
				i = np.mod(((y / self.resolution).decompose().value).astype(np.int32),self.data.shape[0])

			else:

				i = np.mod((y / self.resolution.to(u.rad).value).astype(np.int32),self.data.shape[0])

		else:
			i = None
			j = None
		
		#Call the C backend
		gradient_x,gradient_y = _topology.gradient(self.data,j,i)

		#Return the gradients
		if (x is not None) and (y is not None):

			return gradient_x.reshape(x.shape),gradient_y.reshape(x.shape)

		else:
		
			if save:
				self.gradient_x = gradient_x
				self.gradient_y = gradient_y
		
			return gradient_x,gradient_y

	def hessian(self,x=None,y=None,save=True):
		
		"""
		Computes the hessian of the map and sets the hessian_xx,hessian_yy,hessian_xy attributes accordingly

		:param x: optional, x positions at which to evaluate the hessian
		:type x: array with units

		:param y: optional, y positions at which to evaluate the hessian
		:type y: array with units

		:param save: if True saves the gradient as attrubutes
		:type save: bool.

		:returns: tuple -- (hessian_xx,hessian_yy,hessian_xy)

		>>> test_map = ConvergenceMap.load("map.fit")
		>>> hxx,hyy,hxy = test_map.hessian()

		"""

		if (x is not None) and (y is not None):

			assert x.shape==y.shape,"x and y must have the same shape!"

			#x coordinates
			if type(x)==u.quantity.Quantity:
			
				assert x.unit.physical_type==self.side_angle.unit.physical_type
				j = np.mod(((x / self.resolution).decompose().value).astype(np.int32),self.data.shape[1])

			else:

				j = np.mod((x / self.resolution.to(u.rad).value).astype(np.int32),self.data.shape[1])	

			#y coordinates
			if type(y)==u.quantity.Quantity:
			
				assert y.unit.physical_type==self.side_angle.unit.physical_type
				i = np.mod(((y / self.resolution).decompose().value).astype(np.int32),self.data.shape[0])

			else:

				i = np.mod((y / self.resolution.to(u.rad).value).astype(np.int32),self.data.shape[0])

		else:
			i = None
			j = None

		#Call the C backend
		hessian_xx,hessian_yy,hessian_xy = _topology.hessian(self.data,j,i)
		
		#Return the hessian
		if (x is not None) and (y is not None):

			return hessian_xx.reshape(x.shape),hessian_yy.reshape(x.shape),hessian_xy.reshape(x.shape)

		else:

			if save:
				self.hessian_xx = hessian_xx
				self.hessian_yy = hessian_yy
				self.hessian_xy = hessian_xy

			return hessian_xx,hessian_yy,hessian_xy

	def gradLaplacian(self,x=None,y=None):

		"""
		"""

		if (x is not None) and (y is not None):

			assert x.shape==y.shape,"x and y must have the same shape!"

			#x coordinates
			if type(x)==u.quantity.Quantity:
			
				assert x.unit.physical_type==self.side_angle.unit.physical_type
				j = np.mod(((x / self.resolution).decompose().value).astype(np.int32),self.data.shape[1])

			else:

				j = np.mod((x / self.resolution.to(u.rad).value).astype(np.int32),self.data.shape[1])	

			#y coordinates
			if type(y)==u.quantity.Quantity:
			
				assert y.unit.physical_type==self.side_angle.unit.physical_type
				i = np.mod(((y / self.resolution).decompose().value).astype(np.int32),self.data.shape[0])

			else:

				i = np.mod((y / self.resolution.to(u.rad).value).astype(np.int32),self.data.shape[0])

		else:
			i = None
			j = None

		#Call the C backend
		gl_x,gl_y = _topology.gradLaplacian(self.data,j,i)
		
		#Return the gradient of the laplacian
		if (x is not None) and (y is not None):
			return gl_x.reshape(x.shape),gl_y.reshape(x.shape)
		else:
			return gl_x,gl_y

	################################################################################################################################################

	def pdf(self,thresholds,norm=False):

		"""
		Computes the one point probability distribution function of the convergence map

		:param thresholds: thresholds extremes that define the binning of the pdf
		:type thresholds: array

		:param norm: normalization; if set to a True, interprets the thresholds array as units of sigma (the map standard deviation)
		:type norm: bool.

		:returns: tuple -- (threshold midpoints -- array, pdf normalized at the midpoints -- array)

		:raises: AssertionError if thresholds array is not provided

		>>> test_map = ConvergenceMap.load("map.fit")
		>>> thresholds = np.arange(map.data.min(),map.data.max(),0.05)
		>>> nu,p = test_map.pdf(thresholds)

		"""

		assert thresholds is not None
		midpoints = 0.5 * (thresholds[:-1] + thresholds[1:])

		if norm:
			sigma = self.data.std()
		else:
			sigma = 1.0

		#Compute the histogram
		if self._masked:
			hist,bin_edges = np.histogram(self.data[self._mask],bins=thresholds*sigma,density=True)
		else:
			hist,bin_edges = np.histogram(self.data,bins=thresholds*sigma,density=True)

		#Return
		return midpoints,hist*sigma


	def plotPDF(self,thresholds,norm=False,fig=None,ax=None,**kwargs):

		"""
		Plot the PDF of the map

		"""

		if not matplotlib:
			raise ImportError("matplotlib is not installed, cannot plot the PDF!")

		#Instantiate figure
		if (fig is None) or (ax is None):
			
			self.fig,self.ax = plt.subplots()

		else:

			self.fig = fig
			self.ax = ax

		#Measure the PDF of the pixels
		kappa,pdf = self.pdf(thresholds,norm)

		#Plot the PDF
		self.ax.plot(kappa,pdf,**kwargs)

		#Adjust the labels
		if norm:
			self.ax.set_xlabel(r"$\sigma_{\kappa}$",fontsize=22)
			self.ax.set_ylabel(r"$PDF(\sigma_\kappa)$",fontsize=22)
		else:
			s = self.data.std()
			ax_top = self.ax.twiny()
			ax_top.set_xticks(self.ax.get_xticks())
			ax_top.set_xlim(self.ax.get_xlim())
			ax_top.set_xticklabels([ "{0:.2f}".format(n/s) for n in ax_top.get_xticks() ])

			self.ax.set_xlabel(r"$\kappa$",fontsize=22)
			ax_top.set_xlabel(r"$\kappa/\sigma_\kappa$",fontsize=22)
			self.ax.set_ylabel(r"${\rm PDF}(\kappa)$",fontsize=22)


	################################################################################################################################################


	def peakCount(self,thresholds,norm=False):
		
		"""
		Counts the peaks in the map

		:param thresholds: thresholds extremes that define the binning of the peak histogram
		:type thresholds: array

		:param norm: normalization; if set to a True, interprets the thresholds array as units of sigma (the map standard deviation)
		:type norm: bool.

		:returns: tuple -- (threshold midpoints -- array, differential peak counts at the midpoints -- array)

		:raises: AssertionError if thresholds array is not provided

		>>> test_map = ConvergenceMap.load("map.fit")
		>>> thresholds = np.arange(map.data.min(),map.data.max(),0.05)
		>>> nu,peaks = test_map.peakCount(thresholds)

		"""

		assert thresholds is not None
		midpoints = 0.5 * (thresholds[:-1] + thresholds[1:])

		#Check if the map is masked
		if self._masked:

			if not hasattr(self,"_full_mask"):
				self.maskBoundaries()

			mask_profile = self._full_mask

		else:

			mask_profile = None 

		#Decide if normalizing thresholds
		if norm:
			
			if self._masked:
				sigma = self.data[self._full_mask].std()
			else:
				sigma = self.data.std()

		else:
			sigma = 1.0

		return midpoints,_topology.peakCount(self.data,mask_profile,thresholds,sigma)


	def peakHistogram(self,thresholds,norm=False,fig=None,ax=None,**kwargs):

		"""
		Plot the peak histogram of the map

		"""

		if not matplotlib:
			raise ImportError("matplotlib is not installed, cannot plot the peak histogram!")

		#Instantiate figure
		if (fig is None) or (ax is None):
			
			self.fig,self.ax = plt.subplots()

		else:

			self.fig = fig
			self.ax = ax

		#Count the peaks
		kappa,pk = self.peakCount(thresholds,norm)

		#Plot the peak histogram
		self.ax.plot(kappa,pk*(thresholds[1:]-thresholds[:-1]),**kwargs)
		self.ax.set_yscale("log")

		#Adjust the labels
		if norm:
			self.ax.set_xlabel(r"$\sigma_{\kappa}$",fontsize=22)
			self.ax.set_ylabel(r"$N_{\rm pk}(\sigma_\kappa)$",fontsize=22)
		else:
			s = self.data.std()
			ax_top = self.ax.twiny()
			ax_top.set_xticks(self.ax.get_xticks())
			ax_top.set_xlim(self.ax.get_xlim())
			ax_top.set_xticklabels([ "{0:.2f}".format(n/s) for n in ax_top.get_xticks() ])

			self.ax.set_xlabel(r"$\kappa$",fontsize=22)
			ax_top.set_xlabel(r"$\kappa/\sigma_\kappa$",fontsize=22)
			self.ax.set_ylabel(r"$N_{\rm pk}(\kappa)$",fontsize=22)

	def gaussianPeakHistogram(self,thresholds,norm=False,fig=None,ax=None,**kwargs):
		
		"""
		Plot the Gaussian field approximation (Bond 1987) to the peak histogram

		"""


		if not matplotlib:
			raise ImportError("matplotlib is not installed, cannot plot the peak histogram!")

		#Instantiate figure
		if (fig is None) or (ax is None):
			
			self.fig,self.ax = plt.subplots()

		else:

			self.fig = fig
			self.ax = ax

		#Compute the gradient and laplacian expectation values
		gx,gy = self.gradient()
		hxx,hyy,hxy = self.hessian()

		sigma0 = self.data.std()
		sigma1 = np.sqrt((gx**2+gy**2).mean())
		sigma2 = np.sqrt(((hxx+hyy)**2).mean())

		#Compute special parameters
		g = sigma1**2 / (sigma0*sigma2)
		t = np.sqrt(2)*sigma1/sigma2
		v = 0.5*(thresholds[1:]+thresholds[:-1])
		dv = thresholds[1:] - thresholds[:-1]
		if not(norm):
			v /= sigma0
			dv /= sigma0

		x = v*g

		#Compute G,N
		G = (x**2-g**2)*(1-0.5*sp.erfc(x/np.sqrt(2*(1-g**2)))) + x*(1-g**2)*np.exp(-x**2/(2*(1-g**2)))/np.sqrt(2*np.pi*(1-g**2)) + np.exp(-x**2/(3-2*(g**2)))*(1-0.5*sp.erfc(x/np.sqrt(2*(1-g**2)*(3-2*(g**2)))))/np.sqrt(3-2*(g**2))
		N = np.exp(-0.5*(v**2))*G*dv*(self.data.shape[0]**2)/((t**2)*(2*np.pi)**1.5)

		#Plot the histogram
		self.ax.plot(0.5*(thresholds[1:]+thresholds[:-1]),N,**kwargs)


	def locatePeaks(self,thresholds,norm=False):

		"""
		Locate the peaks in the map

		:param thresholds: thresholds extremes that define the binning of the peak histogram
		:type thresholds: array

		:param norm: normalization; if set to a True, interprets the thresholds array as units of sigma (the map standard deviation)
		:type norm: bool.

		:returns: tuple -- (peak height -- array, peak locations -- array with units)

		:raises: AssertionError if thresholds array is not provided

		>>> test_map = ConvergenceMap.load("map.fit")
		>>> thresholds = np.arange(map.data.min(),map.data.max(),0.05)
		>>> height,positions = test_map.locatePeaks(thresholds)

		"""

		assert thresholds is not None
		if isinstance(thresholds,tuple):
			assert len(thresholds)==2
			thresholds = np.array(thresholds)

		#Check if the map is masked
		if self._masked:

			if not hasattr(self,"_full_mask"):
				self.maskBoundaries()

			mask_profile = self._full_mask

		else:

			mask_profile = None 

		#Decide if normalizing thresholds
		if norm:
			
			if self._masked:
				sigma = self.data[self._full_mask].std()
			else:
				sigma = self.data.std()

		else:
			sigma = 1.0

		#Count the number of pixels between the selected thresholds
		relevant_pixels = ((self.data>=thresholds[0]*sigma) * (self.data<=thresholds[-1]*sigma)).sum()

		#Return the result of the C backend call, scaled to the proper units
		peak_values,peak_locations = _topology.peakLocations(self.data,mask_profile,thresholds,sigma,relevant_pixels)

		return peak_values,(peak_locations*self.resolution).to(self.side_angle.unit)


	def peakDistances(self,thresholds,norm=False):
		
		"""
		Compute the pairwise distance between local maxima on the map

		:param thresholds: thresholds extremes that define the binning of the peak histogram
		:type thresholds: array

		:param norm: normalization; if set to a True, interprets the thresholds array as units of sigma (the map standard deviation)
		:type norm: bool.

		:returns: pairwise distance
		:rtype: quantity

		"""

		#Locate peaks first
		height,loc = self.locatePeaks(thresholds,norm)

		#Measure pairwise distances
		distances = np.sqrt(((loc[None]-loc[:,None])**2).sum(-1))
		i,j = np.indices(distances.shape)

		#Return
		return distances[i>j]


	def peakTwoPCF(self,thresholds,scales,norm=False):

		"""
		Compute the two point function of the peaks on the map

		:param thresholds: thresholds extremes that define the binning of the peak histogram
		:type thresholds: array

		:param scales: radial binning of the 2pcf
		:type scales: quantity

		:param norm: normalization; if set to a True, interprets the thresholds array as units of sigma (the map standard deviation)
		:type norm: bool.

		:returns: (bin centers, peak 2pcf)
		:rtype: tuple

		"""
		raise NotImplementedError


	################################################################################################################################################


	def minkowskiFunctionals(self,thresholds,norm=False):

		"""
		Measures the three Minkowski functionals (area,perimeter and genus characteristic) of the specified map excursion sets

		:param thresholds: thresholds that define the excursion sets to consider
		:type thresholds: array

		:param norm: normalization; if set to a True, interprets the thresholds array as units of sigma (the map standard deviation)
		:type norm: bool.

		:returns: tuple -- (nu -- array, V0 -- array, V1 -- array, V2 -- array) nu are the bins midpoints and V are the Minkowski functionals

		:raises: AssertionError if thresholds array is not provided

		>>> test_map = ConvergenceMap.load("map.fit")
		>>> thresholds = np.arange(-2.0,2.0,0.2)
		>>> nu,V0,V1,V2 = test_map.minkowskiFunctionals(thresholds,norm=True)

		"""

		assert thresholds is not None
		midpoints = 0.5 * (thresholds[:-1] + thresholds[1:])

		#Check if the map is masked
		if self._masked:

			if not hasattr(self,"_full_mask"):
				self.maskBoundaries()

			mask_profile = self._full_mask

		else:

			mask_profile = None 

		#Decide if normalize thresholds or not
		if norm:

			if self._masked:
				sigma = self.data[self._full_mask].std()
			else:
				sigma = self.data.std()
		
		else:
			sigma = 1.0

		#Check if gradient and hessian attributes are available; if not, compute them
		if not (hasattr(self,"gradient_x") and hasattr(self,"gradient_y")):
			self.gradient()

		if not (hasattr(self,"hessian_xx") and hasattr(self,"hessian_yy") and hasattr(self,"hessian_xy")):
			self.hessian()

		#Compute the Minkowski functionals and return them as tuple
		v0,v1,v2 = _topology.minkowski(self.data,mask_profile,self.gradient_x,self.gradient_y,self.hessian_xx,self.hessian_yy,self.hessian_xy,thresholds,sigma)

		return midpoints,v0,v1,v2

	def moments(self,connected=False,dimensionless=False):

		"""
		Measures the first nine moments of the convergence map (two quadratic, three cubic and four quartic)

		:param connected: if set to True returns only the connected part of the moments
		:type connected: bool.

		:param dimensionless: if set to True returns the dimensionless moments, normalized by the appropriate powers of the variance
		:type dimensionless: bool. 

		:returns: array -- (sigma0,sigma1,S0,S1,S2,K0,K1,K2,K3)

		>>> test_map = ConvergenceMap.load("map.fit")
		>>> var0,var1,sk0,sk1,sk2,kur0,kur1,kur2,kur3 = test_map.moments()
		>>> sk0,sk1,sk2 = test_map.moments(dimensionless=True)[2:5]
		>>> kur0,kur1,kur2,kur3 = test_map.moments(connected=True,dimensionless=True)[5:]

		"""

		#First check that the instance has the gradient and hessian attributes; if not, compute them
		if not (hasattr(self,"gradient_x") and hasattr(self,"gradient_y")):
			self.gradient()

		if not (hasattr(self,"hessian_xx") and hasattr(self,"hessian_yy") and hasattr(self,"hessian_xy")):
			self.hessian()

		#Decide if using the full map or only the unmasked region
		if self._masked:

			if not hasattr(self,"_full_mask"):
				self.maskBoundaries()
			
			data = self.data[self._full_mask]
			gradient_x = self.gradient_x[self._full_mask]
			gradient_y = self.gradient_y[self._full_mask]
			hessian_xx = self.hessian_xx[self._full_mask]
			hessian_yy = self.hessian_yy[self._full_mask]
			hessian_xy = self.hessian_xy[self._full_mask]

		else:

			data = self.data
			gradient_x = self.gradient_x
			gradient_y = self.gradient_y
			hessian_xx = self.hessian_xx
			hessian_yy = self.hessian_yy
			hessian_xy = self.hessian_xy

		
		#Quadratic moments
		sigma0 = data.std()
		sigma1 = np.sqrt((gradient_x**2 + gradient_y**2).mean())

		#Cubic moments
		S0 = (data**3).mean()
		S1 = ((data**2)*(hessian_xx + hessian_yy)).mean()
		S2 = ((gradient_x**2 + gradient_y**2)*(hessian_xx + hessian_yy)).mean()

		#Quartic moments
		K0 = (data**4).mean()
		K1 = ((data**3) * (hessian_xx + hessian_yy)).mean()
		K2 = ((data) * (gradient_x**2 + gradient_y**2) * (hessian_xx + hessian_yy)).mean()
		K3 = ((gradient_x**2 + gradient_y**2)**2).mean()

		#Compute connected moments (only quartic affected)
		if connected:
			K0 -= 3 * sigma0**4
			K1 += 3 * sigma0**2 * sigma1**2
			K2 += sigma1**4
			K3 -= 2 * sigma1**4

		
		#Normalize moments to make them dimensionless
		if dimensionless:
			S0 /= sigma0**3
			S1 /= (sigma0 * sigma1**2)
			S2 *= (sigma0 / sigma1**4)

			K0 /= sigma0**4
			K1 /= (sigma0**2 * sigma1**2)
			K2 /= sigma1**4
			K3 /= sigma1**4

			sigma0 /= sigma0
			sigma1 /= sigma1

		#Return the array
		return np.array([sigma0,sigma1,S0,S1,S2,K0,K1,K2,K3])

	################################################################################################################################################

	def powerSpectrum(self,l_edges,scale=None):

		"""
		Measures the power spectrum of the convergence map at the multipole moments specified in the input

		:param l_edges: Multipole bin edges
		:type l_edges: array

		:param scale: scaling to apply to the square of the Fourier pixels before harmonic azimuthal averaging. Must be a function that takes the array of multipole magnitudes as an input and returns an array of real numbers 
		:type scale: callable.

		:returns: (l -- array,Pl -- array) = (binned multipole moments, power spectrum at multipole moments)
		:rtype: tuple.

		>>> test_map = ConvergenceMap.load("map.fit")
		>>> l_edges = np.arange(200.0,5000.0,200.0)
		>>> l,Pl = test_map.powerSpectrum(l_edges)

		"""

		assert not self._masked,"Power spectrum calculation for masked maps is not allowed yet!"
		assert l_edges is not None

		if self.side_angle.unit.physical_type=="length":
			raise NotImplementedError("Power spectrum measurement not implemented yet if side physical unit is length!")

		l = 0.5*(l_edges[:-1] + l_edges[1:])

		#Calculate the Fourier transform of the map with numpy FFT
		ft_map = fftengine.rfft2(self.data)

		#Compute scaling coefficients
		if scale is not None:
			sc = scale(self.getEll())
		else:
			sc = None

		#Compute the power spectrum with the C backend implementation
		power_spectrum = _topology.rfft2_azimuthal(ft_map,ft_map,self.side_angle.to(u.deg).value,l_edges,sc)

		#Output the power spectrum
		return l,power_spectrum

	################################################################################################################################################

	def cross(self,other,statistic="power_spectrum",**kwargs):

		"""
		Measures a cross statistic between maps

		:param other: The other convergence map
		:type other: ConvergenceMap instance

		:param statistic: the cross statistic to measure (default is 'power_spectrum')
		:type statistic: string or callable

		:param kwargs: the keyword arguments are passed to the statistic (when callable)
		:type kwargs: dict.

		:returns: tuple -- (l -- array,Pl -- array) = (multipole moments, cross power spectrum at multipole moments) if the statistic is the power spectrum, otherwise whatever statistic() returns on call

		:raises: AssertionError if the other map has not the same shape as the input one

		>>> test_map = ConvergenceMap.load("map.fit",format=load_fits_default_convergence)
		>>> other_map = ConvergenceMap.load("map2.fit",format=load_fits_default_convergence)
		
		>>> l_edges = np.arange(200.0,5000.0,200.0)
		>>> l,Pl = test_map.cross(other_map,l_edges=l_edges)

		"""
		#The other map must be of the same type as this one
		assert isinstance(other,self.__class__)

		if statistic=="power_spectrum":
			
			assert not (self._masked or other._masked),"Power spectrum calculation for masked maps is not allowed yet!"

			if self.side_angle.unit.physical_type=="length":
				raise NotImplementedError("Power spectrum measurement not implemented yet if side physical unit is length!")

			assert "l_edges" in kwargs
			l_edges = kwargs["l_edges"]

			assert l_edges is not None
			l = 0.5*(l_edges[:-1] + l_edges[1:])

			assert self.side_angle == other.side_angle
			assert self.data.shape == other.data.shape

			#Check if we should scale
			if ("scale" in kwargs) and (kwargs["scale"] is not None):
				sc = kwargs["scale"](self.getEll())
			else:
				sc = None

			#Calculate the Fourier transform of the maps with numpy FFTs
			ft_map1 = fftengine.rfft2(self.data)
			ft_map2 = fftengine.rfft2(other.data)

			#Compute the cross power spectrum with the C backend implementation
			cross_power_spectrum = _topology.rfft2_azimuthal(ft_map1,ft_map2,self.side_angle.to(u.deg).value,l_edges,sc)

			#Output the cross power spectrum
			return l,cross_power_spectrum

		else:

			return statistic(self,other,**kwargs) 


	################################################################################################################################################


	def plotPowerSpectrum(self,l_edges,show_angle=True,angle_unit=u.arcmin,fig=None,ax=None,**kwargs):

		"""
		Plot the power spectrum of the map

		"""

		if not matplotlib:
			raise ImportError("matplotlib is not installed, cannot plot the power spectrum!")

		#Instantiate figure
		if (fig is None) or (ax is None):
			
			self.fig,self.ax = plt.subplots()

		else:

			self.fig = fig
			self.ax = ax

		#Calculate the power spectrum
		l,Pl = self.powerSpectrum(l_edges)

		#Plot
		self.ax.plot(l,l*(l+1.)*Pl/(2*np.pi),**kwargs)
		self.ax.set_xscale("log")
		self.ax.set_yscale("log")

		#Upper scale shows the angle
		if show_angle:
			ax_top = self.ax.twiny()
			ax_top.set_xticks(self.ax.get_xticks())
			ax_top.set_xlim(self.ax.get_xlim())
			ax_top.set_xticklabels(["{0:.2f}".format(((2*np.pi/ell)*u.rad).to(angle_unit).value) for ell in self.ax.get_xticks()])
			ax_top.set_xlabel(r"$\theta$({0})".format(angle_unit.to_string()),fontsize=22)

		#Labels
		self.ax.set_xlabel(r"$\ell$",fontsize=22)
		self.ax.set_ylabel(r"$\ell(\ell+1)P_\ell^{\kappa\kappa}/2\pi$",fontsize=22)


	################################################################################################################################################
	################################################################################################################################################

	def twoPointFunction(self,algorithm="FFT",**kwargs):

		"""
		Computes the two point function of the convergence

		:param algorithm: algorithm used to measure the two point function. Can be "FFT"
		:type algorithm: str.

		:param kwargs: accepted keyword arguments are "theta" to specify the angles at which to measure the 2pcf. If none indicated, the angles are computed automatically
		:type kwargs: dict.

		:returns: (theta,2pcf(theta))
		:rtype: tuple.

		"""

		if algorithm=="FFT":

			#Maximum multipole to use in the calculation
			if "lmax" in kwargs.keys():
				lmax = kwargs["lmax"]
			else:
				lmax = self.lmax

			assert lmax<=self.lmax,"You selected a maximum multipole {0:.3f} which is too high! (Maximum allowed {1:.3f}".format(lmax,self.lmax)

			#Set angles at which to compute the 2pcf
			if "theta" in kwargs.keys():
				theta = kwargs["theta"]
			else:
				theta_min = 2.0*np.pi/lmax
				theta = np.arange(theta_min,self.side_angle.to(u.rad).value,theta_min) * u.rad

			assert theta.unit.physical_type=="angle"

			#Select the multipoles to compute the power spectrum
			l_edges = np.arange(self.lmin,lmax,self.lmin)

			#Compute the power spectrum
			ell,Pell = self.powerSpectrum(l_edges)

			#Use the hankel transform to measure the 2pcf
			two_pcf = fht(0,ell,Pell,theta=theta.to(u.rad).value)[1]

			#Return
			return theta.to(u.arcmin),two_pcf

		else:
			raise NotImplementedError("2PCF algorithm {0} not implemented!".format(algorithm))

	################################################################################################################################################
	################################################################################################################################################

	def countModes(self,l_edges):

		"""
		Counts the available number of modes in multipole space available to each bin specified in l_edges

		:param l_edges: Multipole bin edges
		:type l_edges: array

		:returns: number of modes available to each bin 
		:rtype: array

		:raises: AssertionError if l_edges are not provided

		"""

		assert l_edges is not None

		#Determine the multipole values of each bin in the FFT grid
		ell = self.getEll()

		#Count how many of these pixels fall inside each bin
		modes_on = ell[None] < l_edges[:,None,None]
		modes_ly_0 = modes_on.copy()
		modes_ly_0[:,:,1:] = 0

		#Count the total number of modes, and the number of modes with ly=0 
		num_modes = np.diff(modes_on.sum((1,2)).astype(np.float))
		num_modes_ly_0 = np.diff(modes_ly_0.sum((1,2)).astype(np.float))

		#Return the corrected number of modes that yields the right variance in the Gaussian case
		return num_modes**2/(num_modes+num_modes_ly_0)

	################################################################################################################################################

	def bispectrum(self,l_edges,ratio=0.5,configuration="equilateral",scale=None):

		"""
		Calculates the bispectrum of the map in the equilateral or folded configuration

		:param l_edges: Multipole bin edges: these are the side of the triangle in the equilateral configuration or the base of the triangle in the folded configuration
		:type l_edges: array

		:param ratio: ratio between one of the triangle sides and the base in the folded configuration. Must be between 0 and 1
		:type ratio: float.

		:param configuration: must be either "equilateral" or "folded"
		:type configuration: str.

		:param scale: scaling to apply to the cube of the Fourier pixels before harmonic azimuthal averaging. Must be a function that takes the array of multipole magnitudes as an input and returns an array of real positive numbers
		:type scale: callable.

		:returns: (multipoles, bispectrum at multipoles)
		:rtype: tuple.

		"""

		assert not self._masked,"Bispectrum calculation for masked maps is not allowed yet!"
		assert l_edges is not None

		if self.side_angle.unit.physical_type=="length":
			raise NotImplementedError("Bispectrum measurement not implemented yet if side physical unit is length!")

		#Check folding ratio
		if (configuration=="folded") and not(ratio>0 and ratio<1):
			raise ValueError("Folding ratio should be between 0 and 1!")

		#Multipole edges
		l = 0.5*(l_edges[:-1] + l_edges[1:])

		#Calculate FFT of the map via FFT
		ft_map = fftengine.rfft2(self.data)

		#Scale pixels if scaling is provided
		if scale is not None:
			sc = scale(self.getEll())
			ft_map *= np.cbrt(sc)

		#Calculate bispectrum
		if configuration in ("equilateral","folded"):
			bispectrum = _topology.bispectrum(ft_map,ft_map,ft_map,self.side_angle.to(u.deg).value,l_edges,configuration,ratio)
		else:
			raise NotImplementedError("Bispectrum configuration '{0}' not implemented!".format(configuration))

		#Return
		return l,bispectrum


	################################################################################################################################################

	def smooth(self,scale_angle,kind="gaussian",inplace=False,**kwargs):

		"""
		Performs a smoothing operation on the convergence map

		:param scale_angle: size of the smoothing kernel (must have units)
		:type scale_angle: float.

		:param kind: type of smoothing to be performed. Select "gaussian" for regular Gaussian smoothing in real space or "gaussianFFT" if you want the smoothing to be performed via FFTs (advised for large scale_angle)
		:type kind: str.

		:param inplace: if set to True performs the smoothing in place overwriting the old convergence map
		:type inplace: bool.

		:param kwargs: the keyword arguments are passed to the filter function
		:type kwargs: dict.

		:returns: ConvergenceMap instance (or None if inplace is True)

		"""

		assert not self._masked,"You cannot smooth a masked convergence map!!"

		if kind=="kernelFFT":
			smoothed_data = fftengine.irfft2(scale_angle*fftengine.rfft2(self.data)) 
		else:
			
			assert scale_angle.unit.physical_type==self.side_angle.unit.physical_type

			#Compute the smoothing scale in pixel units
			smoothing_scale_pixel = (scale_angle * self.data.shape[0] / (self.side_angle)).decompose().value

			#Perform the smoothing
			if kind=="gaussian":
				smoothed_data = filters.gaussian_filter(self.data,smoothing_scale_pixel,**kwargs)
		
			elif kind=="gaussianFFT":

				lx = fftengine.fftfreq(self.data.shape[0])
				ly = fftengine.rfftfreq(self.data.shape[1])
				l_squared = lx[:,None]**2 + ly[None,:]**2
				smoothed_data = fftengine.irfft2(np.exp(-0.5*l_squared*(2*np.pi*smoothing_scale_pixel)**2)*fftengine.rfft2(self.data))
		
			else:
				raise NotImplementedError("Smoothing algorithm {0} not implemented!".format(kind))

		#Return the result
		if inplace:
			
			self.data = smoothed_data
			
			#If gradient attributes are present, recompute them
			if (hasattr(self,"gradient_x") or hasattr(self,"gradient_y")):
				self.gradient()

			if (hasattr(self,"hessian_xx") or hasattr(self,"hessian_yy") or hasattr(self,"hessian_xy")):
				self.hessian()
			
			return None

		else:
			
			#Copy the extra attributes as well
			kwargs = dict()
			for attribute in self._extra_attributes:
				kwargs[attribute] = getattr(self,attribute)

			return self.__class__(smoothed_data,self.side_angle,masked=self._masked,**kwargs)


	def __add__(self,rhs):

		"""
		Defines addition operator between ConvergenceMap instances; the convergence values are summed

		:returns: ConvergenceMap instance with the result of the sum

		:raises: AssertionError if the operation cannot be performed

		"""

		if isinstance(rhs,self.__class__):

			assert self.side_angle == rhs.side_angle
			assert self.data.shape == rhs.data.shape

			new_data = self.data + rhs.data

		elif isinstance(rhs,numbers.Number):

			new_data = self.data + rhs

		elif type(rhs) == np.ndarray:

			assert rhs.shape == self.data.shape
			new_data = self.data + rhs

		else:

			raise TypeError("The right hand side cannot be added!!")


		#Copy the extra attributes as well
		kwargs = dict()
		for attribute in self._extra_attributes:
			kwargs[attribute] = getattr(self,attribute)

		return self.__class__(new_data,self.side_angle,masked=self._masked,**kwargs)


	def __mul__(self,rhs):

		"""
		Defines the multiplication operator between ConvergenceMap instances; the convergence values are multiplied (for example the rhs can be a mask...)

		:returns: ConvergenceMap instances with the result of the multiplication

		:raises: AssertionError if the operation cannot be performed

		""" 

		if isinstance(rhs,self.__class__):

			assert self.side_angle == rhs.side_angle
			assert self.data.shape == rhs.data.shape

			new_data = self.data * rhs.data

		elif isinstance(rhs,numbers.Number):

			new_data = self.data * rhs

		elif type(rhs) == np.ndarray:

			assert rhs.shape == self.data.shape
			new_data = self.data * rhs

		else:

			raise TypeError("Cannot multiply by the right hand side!!")

		#Copy the extra attributes as well
		kwargs = dict()
		for attribute in self._extra_attributes:
			kwargs[attribute] = getattr(self,attribute)

		return self.__class__(new_data,self.side_angle,masked=self._masked,**kwargs)


###############################################
##########ConvergenceMap class#################
###############################################

class ConvergenceMap(Spin0):

	"""
	A class that handles 2D convergence maps and allows to compute their topological descriptors (power spectrum, peak counts, minkowski functionals)

	>>> from lenstools import ConvergenceMap 
	>>> from lenstools.defaults import load_fits_default_convergence

	>>> test_map = ConvergenceMap.load("map.fit",format=load_fits_default_convergence)
	>>> imshow(test_map.data)

	"""

#########################################
##########OmegaMap class#################
#########################################

class OmegaMap(Spin0):

	"""
	A class that handles 2D omega maps (curl part of the lensing jacobian) and allows to compute their topological descriptors (power spectrum, peak counts, minkowski functionals)

	"""

###############################################
###################Mask class##################
###############################################

class Mask(ConvergenceMap):

	"""
	A class that handles the convergence masked regions

	"""

	def __init__(self,data,angle,masked=False):

		super(Mask,self).__init__(data,angle,masked)
		assert len(self.data[self.data==0]) + len(self.data[self.data==1]) == reduce(mul,self.data.shape),"The mask must be made of 0 and 1 only!"
		self._masked_fraction = len(self.data[self.data==0]) / reduce(mul,self.data.shape)

	@property
	def maskedFraction(self):

		"""
		Computes the fraction of masked pixels

		"""

		return self._masked_fraction

	def maskBoundaries(self):
		raise AttributeError("Can only call this method on a ConvergenceMap instance!")

	@property
	def boundary(self):
		raise AttributeError("This property is only defined on ConvergenceMap instances!")

######################################################################################################
######################################################################################################

##############
#PhiMap class#
##############

class PhiMap(Spin0):

	"""
	A class that handles 2D lensing potential maps

	"""

#########################
#CMBTemperatureMap class#
#########################

class CMBTemperatureMap(Spin0):

	#Construct sample CMB unlensed map
	@classmethod
	def from_power(cls,angle=3.5*u.deg,npixel=256,seed=None,space="real",powerTT=None,callback="camb_dimensionless",lmax=3500):
		
		"""
		Build an unlensed CMB temperature map from a known TT power spectrum

		:param angle: angular size of the map
		:type angle: quantity

		:param npixel: number of pixels on a side
		:type npixel: int.

		:param seed: random seed for the Fourier coefficients draw
		:type seed: int.

		:param space: must be 'real' or 'fourier'
		:type space: str.

		:param powerTT: name of the file that contains the TT power spectrum. If callback is a callable, powerTT is passed to the callback
		:type powerTT: str.

		:param callback: callback function that computes the TT power spectrum. Can be 'camb_dimensionless' or 'camb_uk' for using camb tabulated power spectra (dimensionless or uK^2 units), None (the identity is used), or callable. If callable, it is called on powerTT and must return (ell,P_TT(ell))
		:type callback: str.

		:returns: temperature map
		:rtype: :py:class:`~lenstools.image.convergence.CMBTemperatureMap`

		"""

		#Random seed
		if seed is not None:
			np.random.seed(seed)

		#CMB lensing routines 
		qlens = Lens()

		#Generate random temperature map
		tfft = qlens.generateTmap(angle,npixel,powerTT,callback,lmax)
		unit = u.uK

		#Return
		if space=="real":
			return cls(fftengine.irfft2(tfft),angle=angle,unit=unit,space="real")
		elif space=="fourier":
			return cls(tfft,angle=angle,unit=unit,space="fourier")
		else:
			raise ValueError("space must be either 'real' or 'fourier'")

	#Add maps taking units into account
	def __add__(self,rhs):

		#Safety check
		assert isinstance(rhs,self.__class__)
		assert self.side_angle==rhs.side_angle

		#Convert units, add
		addedT = self.data + rhs.data*rhs.unit.to(self.unit)

		#Return
		return self.__class__(addedT,self.side_angle,unit=self.unit)

	################################################################################
	################################################################################

	#Switch real/fourier space
	def toReal(self):

		"""
		Switches to real space

		"""

		if self.space=="real":
			pass 
		else:
			self.data = fftengine.irfft2(self.data)
			self.space = "real"

	
	def toFourier(self):

		"""
		Switches to Fourier space via FFT

		"""

		if self.space=="fourier":
			pass 
		else:
			self.data = fftengine.rfft2(self.data)
			self.space="fourier"

	################################################################################

	#Lens the map
	def lens(self,kappa):

		"""
		Lens the CMB temperature map using a kappa map

		:param kappa: convergence map from which the lensing potential is inferred
		:type kappa: :py:class:`~lenstools.image.convergence.ConvergenceMap`

		:returns: lensed temperature map
		:rtype: :py:class:`~lenstools.image.convergence.CMBTemperatureMap`

		"""

		#Safety type check
		assert isinstance(kappa,ConvergenceMap)
		if (self.side_angle!=kappa.side_angle) or (self.data.shape!=kappa.data.shape):
			raise NotImplementedError("Lensing with kappa field of different angle/shape/resolution not implemented yet!")

		#CMB lensing routines 
		qlens = Lens()

		#Lens the map and return
		self.toReal()
		tlens = qlens.lensTmap(self.data,self.side_angle,kappa.data,kappa.side_angle)
		return self.__class__(tlens,self.side_angle,space="real",unit=self.unit)

	################################################################################

	#Add detector effects to the map
	def addDetectorEffects(self,sigmaN=6.0*u.uK*u.arcmin,fwhm=1.4*u.arcmin,seed=None):

		"""
		Add detector effects to the temperature map: convolve with beam, add white noise, deconvolve beam. This operation works in place

		:param sigmaN: rms of the white noise. Must be supplied with temperature x angle units 
		:type sigmaN: quantity

		:param fwhm: full width half maximum of the beam
		:type fwhm: quantity

		:param lmax: maximum multipole to consider
		:type lmax: int.

		:param seed: if not None, seed used to initialize the random number generator
		:type seed: int.

		"""

		#CMB lensing routines
		qlens = Lens()

		#Set random seed
		if seed is not None:
			np.random.seed(seed)

		#Zero out high multipoles (effect of beam + deconvolution on signal)
		tfft = fftengine.rfft2(self.data)
		tfft[qlens._cache["ell"][:,:self.data.shape[0]//2+1]>qlens._cache["lmax"]] = 0.
		self.data = fftengine.irfft2(tfft)

		#rms of the noise in real space
		sigma = (sigmaN / self.resolution).to(self.unit)

		#Generate white noise
		nw = np.random.randn(*self.data.shape)*sigma.value

		#Deconvolve the beam for the noise term
		nw_fft = fftengine.rfft2(nw)
		nw_fft[0,0] = 0.0
		nw_fft *= np.sqrt(qlens._detector(qlens._cache["ell"][:,:self.data.shape[0]//2+1],sigmaN=1.0,fwhm=fwhm.to(u.rad).value,ellmax=qlens._cache["lmax"])) 
		nw = fftengine.irfft2(nw_fft)

		#Add the noise
		self.data += nw

	###########################################
	##########Lensing potential estimation#####
	###########################################

	#Estimate the lensing potential phi using the quadratic estimator on the temperature map
	def estimatePhiFFTQuad(self,powerTT=None,callback="camb_dimensionless",noise_keys=None,lmax=3500,filtering=None):

		"""
		Estimate the Fourier transform of the lensing potential using a temperature quadratic estimator

		:param powerTT: name of the file that contains the lensed theory TT power spectrum. If callback is a callable, powerTT is passed to the callback
		:type powerTT: str.

		:param callback: callback function that computes the TT power spectrum. Can be 'camb_dimensionless' or 'camb_uk' for using camb tabulated power spectra (dimensionless or uK^2 units), None (the identity is used), or callable. If callable, it is called on powerTT and must return (ell,P_TT(ell))
		:type callback: str.

		:param noise_keys: dictionary with noise TT power spectrum specifications
		:type noise_keys: dict.

		:param filtering: filter the map after reconstruction. Can be 'wiener' to apply the wiener filter, or callable. If callable, the function is called on the multipoles and applied to the reconstructed image FFT
		:type filtering: str. or callable

		:returns: FFT of the lensing potential
		:rtype: array

		"""

		#CMB lensing routines 
		qlens = Lens()
		
		#Perform potential estimation with the quadratic estimator (pass the temperature values in uK)
		phifft = qlens.phiTT(fftengine.fft2(self.data)*self.unit.to(u.uK),self.side_angle,powerTT,callback,noise_keys,lmax,filtering)

		#Return
		return phifft

	def estimatePhiQuad(self,powerTT=None,callback="camb_dimensionless",noise_keys=None,lmax=3500,filtering=None):

		"""
		Estimate the lensing potential using a temperature quadratic estimator

		:param powerTT: name of the file that contains the lensed theory TT power spectrum. If callback is a callable, powerTT is passed to the callback
		:type powerTT: str.

		:param callback: callback function that computes the TT power spectrum. Can be 'camb_dimensionless' or 'camb_uk' for using camb tabulated power spectra (dimensionless or uK^2 units), None (the identity is used), or callable. If callable, it is called on powerTT and must return (ell,P_TT(ell))
		:type callback: str.

		:param noise_keys: dictionary with noise TT power spectrum specifications
		:type noise_keys: dict.

		:param filtering: filter the map after reconstruction. Can be 'wiener' to apply the wiener filter, or callable. If callable, the function is called on the multipoles and applied to the reconstructed image FFT
		:type filtering: str. or callable

		:returns: reconstructed lensing potential map
		:rtype: :py:class:`~lenstools.image.convergence.Spin0`

		"""

		#Compute Phi FFT, invert the FFT
		phifft = self.estimatePhiFFTQuad(powerTT,callback,noise_keys,lmax,filtering)
		phi = fftengine.ifft2(phifft)

		#Return
		return PhiMap(phi.real,angle=self.side_angle,unit=u.rad**2)

	#Estimate kappa using the quadratic estimator on the temperature map
	def estimateKappaQuad(self,powerTT=None,callback="camb_dimensionless",noise_keys=None,lmax=3500,filtering=None):

		"""
		Estimate the lensing kappa using a temperature quadratic estimator

		:param powerTT: name of the file that contains the lensed theory TT power spectrum. If callback is a callable, powerTT is passed to the callback
		:type powerTT: str.

		:param callback: callback function that computes the TT power spectrum. Can be 'camb_dimensionless' or 'camb_uk' for using camb tabulated power spectra (dimensionless or uK^2 units), None (the identity is used), or callable. If callable, it is called on powerTT and must return (ell,P_TT(ell))
		:type callback: str.

		:param noise_keys: dictionary with noise TT power spectrum specifications
		:type noise_keys: dict.

		:param filtering: filter the map after reconstruction. Can be 'wiener' to apply the wiener filter, or callable. If callable, the function is called on the multipoles and applied to the reconstructed image FFT
		:type filtering: str. or callable

		:returns: reconstructed convergence map
		:rtype: :py:class:`~lenstools.image.convergence.ConvergenceMap`

		"""

		#CMB lensing routines
		qlens = Lens()

		#Compute Phi FFT, take the laplacian
		phifft = self.estimatePhiFFTQuad(powerTT,callback,noise_keys,lmax,filtering)
		kappafft = phifft*0.5*qlens._cache["ell2"]

		#Invert the FFT
		kappa = fftengine.ifft2(kappafft)

		#Return
		return ConvergenceMap(kappa.real,angle=self.side_angle)

	#Quadratic N0 bias
	def N0Bias(self,l_edges,powerTT=None,callback="camb_dimensionless",noise_keys=None,lmax=3500,output="kappa"):

		"""
		Calculate the quadratic N0 bias on the power spectrum of the kappa/phi estimated with the quadratic TT estimator. The result is referred to the angle and resolution of the image

		:param l_edges: multipole bin edges
		:type l_edges: array.

		:param powerTT: name of the file that contains the lensed theory TT power spectrum. If callback is a callable, powerTT is passed to the callback
		:type powerTT: str.

		:param callback: callback function that computes the TT power spectrum. Can be 'camb_dimensionless' or 'camb_uk' for using camb tabulated power spectra (dimensionless or uK^2 units), None (the identity is used), or callable. If callable, it is called on powerTT and must return (ell,P_TT(ell))
		:type callback: str.

		:param noise_keys: dictionary with noise TT power spectrum specifications
		:type noise_keys: dict.

		:param lmax: maximum multipole to consider
		:type lmax: int.

		:param output: must be 'kappa' or 'phi'
		:type output: str.

		:returns: (bin centers, N0 bias)
		:rtype: tuple.

		"""
		
		#Safety check
		assert output in ("kappa","phi")

		#CMB lensing routines
		qlens = Lens()

		#Perform the estimate
		ell,n0tt = qlens.N0TT(l_edges,self.side_angle,len(self.data),powerTT,callback,noise_keys,lmax)

		#Convert to kappa if requested
		if output=="kappa":
			n0tt *= 0.25*ell**4

		#Return
		return ell,n0tt


