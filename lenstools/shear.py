"""

.. module:: shear 
	:platform: Unix
	:synopsis: This module implements a set of operations which are usually performed on weak lensing shear maps


.. moduleauthor:: Andrea Petri <apetri@phys.columbia.edu>


"""

from __future__ import division

from external import _topology
from topology import ConvergenceMap

import numpy as np
from astropy.io import fits

try:
	import matplotlib.pyplot as plt
	matplotlib = True
except ImportError:
	matplotlib = False

##########################################
#####Default Fits loader##################
##########################################
def load_fits_default(*args):
	"""
	This is the default fits file loader, it assumes that the two components of the shear are stored in two different image FITS files, which have an ANGLE keyword in the header

	:param gamma_file: Name of the FITS file that contains the shear map
	:type gamma1_file: str.

	:returns: tuple -- (angle,ndarray -- gamma; gamma[0] is the gamma1 map, gamma[1] is the gamma2 map); the maps must follow matrix ordering, i.e. the first axis (0) is y and the second axis (1) is x. This matters for the E/B mode decomposition 

	:raises: IOError if the FITS files cannot be opened or do not exist

	"""

	#Open the files
	gamma_file = fits.open(args[0])

	#Read the ANGLE keyword from the header
	angle = gamma_file[0].header["ANGLE"]

	#Create the array with the shear map
	gamma = gamma_file[0].data.astype(np.float)

	#Close files and return
	gamma_file.close()

	return angle,gamma


##########################################
########ShearMap class####################
##########################################

class ShearMap(object):

	"""
	A class that handles 2D shear maps and allows to perform a set of operations on them

	>>> from lenstools.shear import ShearMap
	
	>>> test = ShearMap.fromfilename("shear.fit",loader=load_fits_default)
	>>>test.side_angle
	1.95
	>>>test.gamma
	#The actual map values

	"""

	def __init__(self,gamma,angle):

		self.gamma = gamma
		self.side_angle = angle

	@classmethod
	def fromfilename(cls,*args,**kwargs):
		
		"""
		This class method allows to read the map from a data file; the details of the loading are performed by the loader function. The only restriction to this function is that it must return a tuple (angle,gamma)

		:param args: The positional arguments that are to be passed to the loader (typically the file name)

		:param kwargs: Only one keyword is accepted "loader" is a pointer to the previously defined loader method (the default is load_fits_default above)
		
		"""

		if not("loader" in kwargs.keys()):
			loader = load_fits_default
		else:
			loader = kwargs["loader"]

		angle,gamma = loader(*args)
		return cls(gamma,angle)

	def sticks(self,ax,pixel_step=10):

		"""
		Draw the ellipticity map using the shear components

		:param ax: ax on which to draw the ellipticity field
		:type ax: matplotlib ax object

		:param pixel_step: One arrow will be drawn every pixel_step pixels to avoid arrow overplotting
		:type pixel_step: int.

		>>> import matplotlib.pyplot as plt
		>>> test = ShearMap.fromfilename("shear.fit",loader=load_fits_default)
		>>> fig,ax = plt.subplots()
		>>> test.sticks(ax,pixel_step=50)

		"""

		if not(matplotlib):
			print("matplotlib is not installed, please install it!!")
			return None

		x,y = np.meshgrid(np.arange(0,self.gamma.shape[2],pixel_step),np.arange(0,self.gamma.shape[1],pixel_step))

		#Translate shear components into sines and cosines
		cos_2_phi = self.gamma[0] / np.sqrt(self.gamma[0]**2 + self.gamma[1]**2)
		sin_2_phi = self.gamma[1] / np.sqrt(self.gamma[0]**2 + self.gamma[1]**2)

		#Compute stick directions
		cos_phi = np.sqrt(0.5*(1.0 + cos_2_phi)) * np.sign(sin_2_phi)
		sin_phi = np.sqrt(0.5*(1.0 - cos_2_phi))

		#Draw map using matplotlib quiver
		if ax is None:
			fig,ax = plt.subplots()

		ax.quiver(x*self.side_angle/self.gamma.shape[2],y*self.side_angle/self.gamma.shape[1],cos_phi[x,y],sin_phi[x,y],headwidth=0,units="height",scale=x.shape[0])

		return ax




	def decompose(self,l_edges,keep_fourier=False):

		"""
		Decomposes the shear map into its E and B modes components and returns the respective power spectal densities at the specified multipole moments

		:param l_edges: Multipole bin edges
		:type l_edges: array

		:param keep_fourier: If set to True, holds the Fourier transforms of the E and B mode maps into the E and B attributes of the ShearMap instance
		:type keep_fourier: bool. 

		:returns: :returns: tuple -- (l -- array,P_EE,P_BB,P_EB -- arrays) = (multipole moments, EE,BB power spectra and EB cross power)

		>>> test_map = ShearMap.fromfilename("shear.fit",loader=load_fits_default)
		>>> l_edges = np.arange(300.0,5000.0,200.0)
		>>> l,EE,BB,EB = test_map.decompose(l_edges)

		"""

		#Perform Fourier transforms
		ft_gamma1 = np.fft.rfft2(self.gamma[0])
		ft_gamma2 = np.fft.rfft2(self.gamma[1])

		#Compute frequencies
		lx = np.fft.rfftfreq(ft_gamma1.shape[0])
		ly = np.fft.fftfreq(ft_gamma1.shape[0])

		#Safety check
		assert len(lx)==ft_gamma1.shape[1]
		assert len(ly)==ft_gamma1.shape[0]

		#Compute sines and cosines of rotation angles
		l_squared = lx[np.newaxis,:]**2 + ly[:,np.newaxis]**2
		l_squared[0,0] = 1.0

		sin_2_phi = 2.0 * lx[np.newaxis,:] * ly[:,np.newaxis] / l_squared
		cos_2_phi = (lx[np.newaxis,:]**2 - ly[:,np.newaxis]**2) / l_squared

		#Compute E and B components
		ft_E = cos_2_phi * ft_gamma1 + sin_2_phi * ft_gamma2
		ft_B = -1.0 * sin_2_phi * ft_gamma1 + cos_2_phi * ft_gamma2

		ft_E[0,0] = 0.0
		ft_B[0,0] = 0.0

		assert ft_E.shape == ft_B.shape
		assert ft_E.shape == ft_gamma1.shape

		#Compute and return power spectra
		l = 0.5*(l_edges[:-1] + l_edges[1:])
		P_EE = _topology.rfft2_azimuthal(ft_E,ft_E,self.side_angle,l_edges)
		P_BB = _topology.rfft2_azimuthal(ft_B,ft_B,self.side_angle,l_edges)
		P_EB = _topology.rfft2_azimuthal(ft_E,ft_B,self.side_angle,l_edges)

		if keep_fourier:
			self.fourier_E = ft_E
			self.fourier_B = ft_B

		return l,P_EE,P_BB,P_EB



	def convergence(self):
		
		"""
		Reconstructs the convergence from the E component of the shear

		:returns: new ConvergenceMap instance 

		"""

		#Compute Fourier transforms if it wasn't done before
		if not hasattr(self,"fourier_E"):
			l_edges = np.array([200.0,400.0])
			l,EE,BB,EB = self.decompose(l_edges,keep_fourier=True)

		#Invert the Fourier transform to go back to real space
		conv = np.fft.irfft2(self.fourier_E)

		#Return the ConvergenceMap instance
		return ConvergenceMap(conv,self.side_angle)



