from ..convergence import Spin0,ConvergenceMap
from ..shear import Spin1,Spin2,ShearMap

import time
import logging
import re

from operator import mul
from functools import reduce

import numpy as np
from scipy.spatial import cKDTree as KDTree

try:
	import matplotlib.pyplot as plt
	matplotlib = plt
except ImportError:
	matplotlib = None

#FFT engines
from numpy.fft import fftfreq,rfft2,irfft2

try:
	from numpy.fft import rfftfreq
except ImportError:
	from ..utils import rfftfreq

from astropy.cosmology import w0waCDM
import astropy.units 
from astropy.units import km,s,Mpc,rad,deg,dimensionless_unscaled,quantity

#Try to import the FITSIO library for optimal FITS images reading
try:
	from fitsio import FITS as fitsio
	fitsio = fitsio
	from astropy.io import fits

except ImportError:
	from astropy.io import fits
	fitsio = None


###########################################################
#################Plane class###############################
###########################################################

class Plane(Spin0):


	def __init__(self,data,angle,redshift=2.0,cosmology=None,comoving_distance=None,unit=rad**2,num_particles=None,masked=False):

		#Sanity check
		assert (cosmology is not None) or (comoving_distance is not None),"cosmology and comoving_distance cannot be both None!!"

		super(Plane,self).__init__(data,angle,masked=masked,redshift=redshift,cosmology=cosmology,comoving_distance=comoving_distance,unit=unit,num_particles=num_particles)
		self.redshift = redshift
		self.cosmology = cosmology
		self.unit = unit

		if num_particles is not None:
			self.num_particles = num_particles
		else:
			self.num_particles = -1

		#If a comoving distance is provided, use that; otherwise it needs to be computed from the astropy cosmology instance
		if comoving_distance is not None:		
			
			assert comoving_distance.unit.physical_type=="length"
			self.comoving_distance = comoving_distance
		
		else:
			self.comoving_distance = cosmology.comoving_distance(redshift)

		if data.dtype in [np.float,np.float32]:
			self.space = "real"
		elif data.dtype==np.complex:
			self.space = "fourier"
		else:
			raise TypeError("data type not supported!")


	def angular(self,length_scale):

		"""
		Given a comoving spatial length scale, this method computes the corresponding angular length scale on the plane in both physical and pixel units

		:param length_scale: comoving length scale to be converted into angular units
		:type length_scale: float with units

		:returns: tuple(length_scale in degrees,length_scale in pixels)

		"""

		#Length scale must have units of length
		assert length_scale.unit.physical_type=="length"
		assert self.side_angle.unit.physical_type=="angle"

		angle_scale = (length_scale/self.comoving_distance).decompose().value * rad
		pixel_scale = (angle_scale / self.resolution).decompose().value

		#return
		return angle_scale.to(deg),pixel_scale


	def save(self,filename,format=None,double_precision=False):

		"""
		Saves the Plane to an external file, of which the format can be specified (only fits implemented so far)

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

			#A cosmology instance should be available in order to save in FITS format
			assert self.cosmology is not None
		
			#Create the hdu
			if self.space=="real":
				
				if double_precision:
					hdu = fits.PrimaryHDU(self.data)
				else:
					hdu = fits.PrimaryHDU(self.data.astype(np.float32))
			
			elif self.space=="fourier":
				
				hdu = fits.PrimaryHDU(self.data.real)
				hdu1 = fits.ImageHDU(self.data.imag)
			
			else:
				raise TypeError("data type not supported!")


			#Generate a header
			hdu.header["H0"] = (self.cosmology.H0.to(km/(s*Mpc)).value,"Hubble constant in km/s*Mpc")
			hdu.header["h"] = (self.cosmology.h,"Dimensionless Hubble constant")
			hdu.header["OMEGA_M"] = (self.cosmology.Om0,"Dark Matter density")
			hdu.header["OMEGA_L"] = (self.cosmology.Ode0,"Dark Energy density")
			hdu.header["W0"] = (self.cosmology.w0,"Dark Energy equation of state")
			hdu.header["WA"] = (self.cosmology.wa,"Dark Energy running equation of state")

			hdu.header["Z"] = (self.redshift,"Redshift of the lens plane")
			hdu.header["CHI"] = (hdu.header["h"] * self.comoving_distance.to(Mpc).value,"Comoving distance in Mpc/h")

			if self.side_angle.unit.physical_type=="angle":
				hdu.header["ANGLE"] = (self.side_angle.to(deg).value,"Side angle in degrees")
			elif self.side_angle.unit.physical_type=="length":
				hdu.header["SIDE"] = (self.side_angle.to(Mpc).value*self.cosmology.h,"Side length in Mpc/h")

			hdu.header["NPART"] = (float(self.num_particles),"Number of particles on the plane")
			hdu.header["UNIT"] = (self.unit.to_string(),"Measure units of the pixel values") 

			#Save the plane
			if self.space=="real":
				hdulist = fits.HDUList([hdu])
			else:
				hdulist = fits.HDUList([hdu,hdu1])

			hdulist.writeto(filename,clobber=True)

		else:
			raise ValueError("Format {0} not implemented yet!!".format(format))


	@classmethod
	def load(cls,filename,format=None,init_cosmology=True):

		"""
		Loads the Plane from an external file, of which the format can be specified (only fits implemented so far)

		:param filename: name of the file from which to load the plane
		:type filename: str.

		:param format: format of the file, only FITS implemented so far; if None, it's detected automatically from the filename
		:type format: str.

		:param init_cosmology: if True, instantiates the cosmology attribute of the PotentialPlane
		:type init_cosmology: bool.

		:returns: PotentialPlane instance that wraps the data contained in the file

		"""

		if format is None:
			
			extension = filename.split(".")[-1]
			if extension in ["fit","fits"]:
				format="fits"
			else:
				raise IOError("File format not recognized from extension '{0}', please specify it manually".format(extension))


		if format=="fits":

			#Read the FITS file with the plane information (if there are two HDU's the second one is the imaginary part)
			if fitsio is not None:
				hdu = fitsio(filename)
			else:
				hdu = fits.open(filename)
			
			if len(hdu)>2:
				raise ValueError("There are more than 2 HDUs, file format unknown")

			if fitsio is not None:
				header = hdu[0].read_header()
			else:
				header = hdu[0].header

			#Retrieve the info from the header (handle old FITS header format too)
			try:
				hubble = header["H0"] * (km/(s*Mpc))
				h = header["h"]
			except:
				hubble = header["H_0"] * (km/(s*Mpc))
				h = hubble.value / 100

			Om0 = header["OMEGA_M"]
			Ode0 = header["OMEGA_L"]

			try:
				w0 = header["W0"]
				wa = header["WA"]
			except:
				w0 = header["W_0"]
				wa = header["W_A"]
			
			redshift = header["Z"]
			comoving_distance = (header["CHI"] / h) * Mpc

			if "SIDE" in header.keys():
				angle = header["SIDE"] * Mpc / h
			elif "ANGLE" in header.keys():
				angle = header["ANGLE"] * deg
			else:
				angle = ((header["RES_X"] * header["NAXIS1"] / header["CHI"]) * rad).to(deg)

			#Build the cosmology object if options directs
			if init_cosmology:
				cosmology = w0waCDM(H0=hubble,Om0=Om0,Ode0=Ode0,w0=w0,wa=wa)
			else:
				cosmology = None

			#Read the number of particles, if present
			try:
				num_particles = header["NPART"]
			except:
				num_particles = None

			#Read the units if present
			try:
				unit_string = header["UNIT"]
				name,exponent = re.match(r"([a-zA-Z]+)([0-9])?",unit_string).groups()
				unit = getattr(astropy.units,name)
				if exponent is not None:
					unit **= exponent
			except AttributeError:
				unit = dimensionless_unscaled
			except (ValueError,KeyError):
				unit = rad**2

			#Instantiate the new PotentialPlane instance
			if fitsio is not None:

				if len(hdu)==1:
					new_plane = cls(hdu[0].read(),angle=angle,redshift=redshift,comoving_distance=comoving_distance,cosmology=cosmology,unit=unit,num_particles=num_particles)
				else:
					new_plane = cls(hdu[1].read() + 1.0j*hdu[1].read(),angle=angle,redshift=redshift,comoving_distance=comoving_distance,cosmology=cosmology,unit=unit,num_particles=num_particles)

			else:
			
				if len(hdu)==1:
					new_plane = cls(hdu[0].data.astype(np.float64),angle=angle,redshift=redshift,comoving_distance=comoving_distance,cosmology=cosmology,unit=unit,num_particles=num_particles)
				else:
					new_plane = cls((hdu[0].data + 1.0j*hdu[1].data).astype(np.complex128),angle=angle,redshift=redshift,comoving_distance=comoving_distance,cosmology=cosmology,unit=unit,num_particles=num_particles)

			#Close the FITS file and return
			hdu.close()
			return new_plane

		else:
			raise ValueError("Format {0} not implemented yet!!".format(format))


	def randomRoll(self,seed=None,lmesh=None):

		"""
		Randomly shifts the plane along its axes, enforcing periodic boundary conditions

		:param seed: random seed with which to initialize the generator
		:type seed: int.

		:param lmesh: the FFT frequency meshgrid (lx,ly) necessary for the calculations in fourier space; if None, a new one is computed from scratch (must have the appropriate dimensions)
		:type lmesh: array

		"""

		now = time.time()
		last_timestamp = now

		if seed is not None:
			np.random.seed(seed)

		if self.space=="real":

			#Roll in real space
			self.data = np.roll(np.roll(self.data,np.random.randint(0,self.data.shape[0]),axis=0),np.random.randint(0,self.data.shape[1]),axis=1)	
		
		elif self.space=="fourier":

			#Rolling in Fourier space is just multiplying by phases
			if lmesh is None:
				l = np.array(np.meshgrid(rfftfreq(self.data.shape[0]),fftfreq(self.data.shape[0])))
			else:
				l = lmesh

			#Timestamp
			now = time.time()
			logging.debug("l meshgrid initialized in {0:.3f}s".format(now-last_timestamp))
			last_timestamp = now 

			random_shift = np.random.randint(0,self.data.shape[0],size=2)
			self.data *= np.exp(2.0j*np.pi*np.tensordot(random_shift,l,axes=(0,0)))

			#Timestamp
			now = time.time()
			logging.debug("Phase multiplication completed in {0:.3f}s".format(now-last_timestamp))
			last_timestamp = now 


		else:
			raise ValueError("space must be either real or fourier!")


	def toReal(self):

		"""
		Switches to real space

		"""

		assert self.space=="fourier","We are already in real space!!"
		self.data = irfft2(self.data)
		self.space = "real"

	
	def toFourier(self):

		"""
		Switches to Fourier space

		"""

		assert self.space=="real","We are already in fourier space!!"
		self.data = rfft2(self.data)
		self.space="fourier"


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
		if type(x)==quantity.Quantity:
			
			assert x.unit.physical_type=="angle"

			#Check if the resolution units are length units
			if self.resolution.unit.physical_type=="length":
				x = x.to(rad).value*self.comoving_distance 

			j = np.mod(((x / self.resolution).decompose().value).astype(np.int32),self.data.shape[1])

		else:

			j = np.mod((x / self.resolution.to(rad).value).astype(np.int32),self.data.shape[1])	

		#y coordinates
		if type(y)==quantity.Quantity:
			
			assert y.unit.physical_type=="angle"

			#Check if the resolution units are length units
			if self.resolution.unit.physical_type=="length":
				y = y.to(rad).value*self.comoving_distance

			i = np.mod(((y / self.resolution).decompose().value).astype(np.int32),self.data.shape[0])

		else:

			i = np.mod((y / self.resolution.to(rad).value).astype(np.int32),self.data.shape[0])

		#Return the map values at the specified coordinates
		return self.data[i,j]


###########################################################
#################DensityPlane class########################
###########################################################

class DensityPlane(Plane):

	"""
	Class handler of a lens density plane, inherits from the parent Plane class; additionally it defines redshift and comoving distance attributes that are needed for ray tracing operations

	"""

	def potential(self,lmesh=None):

		"""
		Computes the lensing potential from the density plane solving the Poisson equation via FFTs

		:param lmesh: the FFT frequency meshgrid (lx,ly) necessary for the calculations in fourier space; if None, a new one is computed from scratch (must have the appropriate dimensions)
		:type lmesh: array

		:returns: PotentialPlane instance with the computed lensing potential

		"""

		if self.side_angle.unit.physical_type=="length":
			raise NotImplementedError("potential calculations for physical planes not implemented yet")

		#Initialize l meshgrid
		if lmesh is None:
			l = np.array(np.meshgrid(rfftfreq(self.data.shape[0]),fftfreq(self.data.shape[0])))
		else:
			l = lmesh

		#Compute the magnitude squared of the wavenumber
		l_squared = l[0]**2 + l[1]**2
		l_squared[0,0] = 1.0

		#Go with the FFTs
		if self.space=="real":
			density_ft = rfft2(self.data)
		elif self.space=="fourier":
			density_ft = self.data.copy()
		else:
			raise ValueError("space must be either real or fourier!")


		#Invert the laplacian
		density_ft *= -2.0*((self.resolution.to(rad).value)**2) / (l_squared * ((2.0*np.pi)**2))
		density_ft[0,0] = 0.0

		#Instantiate the new PotentialPlane
		return PotentialPlane(data=irfft2(density_ft),angle=self.side_angle,redshift=self.redshift,comoving_distance=self.comoving_distance,cosmology=self.cosmology,num_particles=self.num_particles,unit=rad**2)


###########################################################
###############PotentialPlane class########################
###########################################################

class PotentialPlane(Plane):

	"""
	Class handler of a lens potential plane, inherits from the parent Plane class; additionally it defines redshift and comoving distance attributes that are needed for ray tracing operations

	"""

	
	def deflectionAngles(self,x=None,y=None,lmesh=None):

		"""
		Computes the deflection angles for the given lensing potential by taking the gradient of the potential map; it is also possible to proceed with FFTs

		:param x: optional; if not None, compute the deflection angles only for rays hitting the lens at the particular x positions (mainly for speedup in case there are less light rays than the plane resolution allows; must proceed in real space to allow speedup)
		:type x: array with units

		:param y: optional; if not None, compute the deflection angles only for rays hitting the lens at the particular y positions (mainly for speedup in case there are less light rays than the plane resolution allows; must proceed in real space to allow speedup)
		:type y: array with units

		:param lmesh: the FFT frequency meshgrid (lx,ly) necessary for the calculations in fourier space; if None, a new one is computed from scratch (must have the appropriate dimensions)
		:type lmesh: array

		:returns: DeflectionPlane instance, or array with deflections of rays hitting the lens at (x,y)

		"""

		now = time.time()
		last_timestamp = now

		if self.space=="real":

			#Scale x and y to lengths in case this is a physical plane
			if self.side_angle.unit.physical_type=="length" and (x is not None) and (y is not None):
				x = x.to(rad).value * self.comoving_distance
				y = y.to(rad).value * self.comoving_distance
			
			#Compute the gradient of the potential map
			deflection_x,deflection_y = self.gradient(x,y)
			deflection = np.array([deflection_x,deflection_y])
		
		elif self.space=="fourier":

			#It doesn't make sense to select (x,y) when we proceed with FFTs
			assert (x is None) and (y is None),"It doesn't make sense to select (x,y) when we proceed with FFTs!"

			#Compute deflections in fourier space
			if lmesh is None:
				l = np.array(np.meshgrid(rfftfreq(self.data.shape[0]),fftfreq(self.data.shape[0])))
			else:
				l = lmesh

			#Timestamp
			now = time.time()
			logging.debug("l meshgrid initialized in {0:.3f}s".format(now-last_timestamp))
			last_timestamp = now 

			ft_deflection = 2.0*np.pi*1.0j * l * self.data

			#Timestamp
			now = time.time()
			logging.debug("Phase multiplications completed in {0:.3f}s".format(now-last_timestamp))
			last_timestamp = now 

			#Go back in real space
			deflection = irfft2(ft_deflection)

			#Timestamp
			now = time.time()
			logging.debug("Inverse FFTs completed in {0:.3f}s".format(now-last_timestamp))
			last_timestamp = now 

		else:
			raise ValueError("space must be either real or fourier!")

		#Scale to units
		deflection *= self.unit
		deflection /= self.resolution
		
		if self.side_angle.unit.physical_type=="length":
			deflection *= self.comoving_distance
			deflection /= rad
		
		assert deflection.unit.physical_type=="angle"
		deflection = deflection.to(rad)

		if (x is not None) and (y is not None):
			
			#If x and y are specified, return the deflections only at those particular points
			return deflection

		else:
		
			#Otherwise return the whole DeflectionPlane instance if we computed the entire mesh
			return DeflectionPlane(deflection,angle=self.side_angle,redshift=self.redshift,comoving_distance=self.comoving_distance,cosmology=self.cosmology,unit=rad)


	def shearMatrix(self,x=None,y=None,lmesh=None):

		"""
		Computes the shear matrix for the given lensing potential; it is also possible to proceed with FFTs

		:param x: optional; if not None, compute the shear matrix only for rays hitting the lens at the particular x positions (mainly for speedup in case there are less light rays than the plane resolution allows; must proceed in real space to allow speedup)
		:type x: array with units

		:param y: optional; if not None, compute the shear matrix only for rays hitting the lens at the particular y positions (mainly for speedup in case there are less light rays than the plane resolution allows; must proceed in real space to allow speedup)
		:type y: array with units

		:param lmesh: the FFT frequency meshgrid (lx,ly) necessary for the calculations in fourier space; if None, a new one is computed from scratch (must have the appropriate dimensions)
		:type lmesh: array

		:returns: ShearTensorPlane instance, or array with deflections of rays hitting the lens at (x,y)

		"""

		now = time.time()
		last_timestamp = now

		if self.space=="real":

			#Scale x and y to lengths in case this is a physical plane
			if self.side_angle.unit.physical_type=="length" and (x is not None) and (y is not None):
				x = x.to(rad).value * self.comoving_distance
				y = y.to(rad).value * self.comoving_distance
			
			#Compute the second derivatives
			tensor = np.array(self.hessian(x,y))

		elif self.space=="fourier":

			#It doesn't make sense to select (x,y) when we proceed with FFTs
			assert (x is None) and (y is None),"It doesn't make sense to select (x,y) when we proceed with FFTs!"

			#Compute deflections in fourier space
			if lmesh is None:
				lx,ly = np.array(np.meshgrid(rfftfreq(self.data.shape[0]),fftfreq(self.data.shape[0])))
			else:
				lx,ly = lmesh

			#Compute deflections in fourier space
			ft_tensor_xx = -(2.0*np.pi)**2 * self.data * (lx**2)
			ft_tensor_xy = -(2.0*np.pi)**2 * self.data * (lx*ly)
			ft_tensor_yy = -(2.0*np.pi)**2 * self.data * (ly**2)

			#Go back in real space
			tensor_xx = irfft2(ft_tensor_xx)
			tensor_xy = irfft2(ft_tensor_xy)
			tensor_yy = irfft2(ft_tensor_yy)

			tensor = np.array([tensor_xx,tensor_yy,tensor_xy])

		else:
			raise ValueError("space must be either real or fourier!")


		#Scale units
		tensor *= self.unit
		tensor /= self.resolution**2

		if self.side_angle.unit.physical_type=="length":
			tensor *= self.comoving_distance**2
			tensor /= rad**2

		assert tensor.unit.physical_type=="dimensionless"
		tensor = tensor.decompose().value

		#Return the ShearTensorPlane instance
		if (x is not None) and (y is not None):

			return tensor

		else:
			return ShearTensorPlane(tensor,angle=self.side_angle,redshift=self.redshift,comoving_distance=self.comoving_distance,cosmology=self.cosmology,unit=dimensionless_unscaled)


	def density(self,x=None,y=None):

		"""
		Computes the projected density fluctuation by taking the laplacian of the potential; useful to check if the potential is reasonable

		:param x: optional; if not None, compute the density only for rays hitting the lens at the particular x positions (mainly for speedup in case there are less light rays than the plane resolution allows; must proceed in real space to allow speedup)
		:type x: array with units

		:param y: optional; if not None, compute the density only for rays hitting the lens at the particular y positions (mainly for speedup in case there are less light rays than the plane resolution allows; must proceed in real space to allow speedup)
		:type y: array with units

		:returns: DensityPlane instance with the density fluctuation data (if x and y are None), or numpy array with the same shape as x and y 

		"""

		#Compute the laplacian
		if self.space=="real":

			#Scale x and y to lengths in case this is a physical plane
			if self.side_angle.unit.physical_type=="length" and (x is not None) and (y is not None):
				x = x.to(rad).value * self.comoving_distance
				y = y.to(rad).value * self.comoving_distance			
			
			logging.debug("Computing hessian...")
			hessian_xx,hessian_yy,hessian_xy = self.hessian(x,y)
			
			logging.debug("Computing laplacian...")
			laplacian = hessian_xx + hessian_yy
			
			logging.debug("Laplacian calculation completed")

		elif self.space=="fourier":

			ly,lx = np.meshgrid(fftfreq(self.data.shape[0]),rfftfreq(self.data.shape[0]),indexing="ij")
			ft_laplacian = -1.0 * (2.0*np.pi)**2 * (lx**2 + ly**2) * self.data
			laplacian = irfft2(ft_laplacian) 			

		else:
			raise ValueError("space must be either real or fourier!")

		#Scale the units
		laplacian *= self.unit
		laplacian /= self.resolution**2

		if self.side_angle.unit.physical_type=="length":
			laplacian *= self.comoving_distance**2
			laplacian /= rad**2

		assert laplacian.unit.physical_type=="dimensionless"
		laplacian = laplacian.decompose().value

		#The density is twice the trace of the hessian
		if (x is not None) and (y is not None):
			return 0.5*laplacian
		else:
			return DensityPlane(0.5*laplacian,angle=self.side_angle,cosmology=self.cosmology,redshift=self.redshift,comoving_distance=self.comoving_distance,num_particles=self.num_particles,unit=dimensionless_unscaled)


#############################################################
################DeflectionPlane class########################
#############################################################

class DeflectionPlane(Spin1):

	"""
	Class handler of a lens deflection plane, inherits from the parent Spin1 class and holds the values of the deflection angles of the light rays that cross a potential plane

	"""

	def __init__(self,data,angle,redshift=2.0,comoving_distance=None,cosmology=None,unit=rad):

		#Sanity check
		assert (cosmology is not None) or (comoving_distance is not None),"cosmology and comoving_distance cannot be both none!!"

		super(self.__class__,self).__init__(data,angle)
		self.redshift = redshift
		self.unit = unit

		#If a comoving distance is provided, use that; otherwise it needs to be computed from the astropy cosmology instance
		if comoving_distance is not None:		
			
			assert comoving_distance.unit.physical_type=="length"
			self.comoving_distance = comoving_distance
		
		else:
			self.comoving_distance = cosmology.comoving_distance(redshift)


	def jacobian(self):

		"""
		Computes the jacobian of the deflection angles, useful to compute shear and convergence; units are handled properly

		:returns: the jacobian of the deflection field in array form, of shape (4,:,:) where the four components are, respectively, xx,xy,yx,yy

		"""

		jac = self.gradient() * self.unit / self.resolution
		
		if self.side_angle.unit.physical_type=="length":
			jac *= self.comoving_distance
			jac /= rad

		assert jac.unit.physical_type=="dimensionless"
		return jac.decompose().value


	def convergence(self,precision="first"):

		"""
		Computes the convergence from the deflection angles by taking the appropriate components of the jacobian

		:param precision: if "first" computes the convergence at first order in the lensing potential (only one implemented so far)
		:type precision: str.

		:returns: ConvergenceMap instance with the computed convergence values

		"""

		#Compute the jacobian and the convergence by tracing it
		jacobian = self.jacobian()
		convergence = -0.5*(jacobian[0]+jacobian[3])

		#Instantiate the convergence map
		return ConvergenceMap(convergence,angle=self.side_angle)


	def omega(self,precision="first"):

		"""
		Computes the omega field (i.e. the real space B mode) from the deflection angles by taking the appropriate components of the jacobian

		:param precision: if "first" computes omega at first order in the lensing potential (only one implemented so far)
		:type precision: str.

		:returns: Spin0 instance with the computed omega values

		"""

		#Compute the jacobian and the convergence by tracing it
		jacobian = self.jacobian()
		omega = -0.5*(jacobian[2]-jacobian[1])

		#Instantiate the convergence map
		return Spin0(omega,angle=self.side_angle)


	def shear(self,precision="first"):

		"""
		Computes the shear from the deflection angles by taking the appropriate components of the jacobian

		:param precision: if "first" computes the shear at first order in the lensing potential (only one implemented so far)
		:type precision: str. 

		:returns: ShearMap instance with the computed convergence values

		"""

		#Compute the jacobian and the shear
		jacobian = self.jacobian()
		shear = np.array([0.5*(jacobian[3]-jacobian[0]),-0.5*(jacobian[1]+jacobian[2])])

		#Instantiate the shear map
		return ShearMap(shear,angle=self.side_angle)


############################################################
#############ShearTensorPlane class#########################
############################################################

class ShearTensorPlane(Spin2):

	"""
	Class handler of a plane of shear matrices, inherits from the parent Spin2 class and holds the 3 values of the symmetric shear matrix (2 diagonal + 1 off diagonal), for each pixel

	"""

	def __init__(self,data,angle,redshift=2.0,comoving_distance=None,cosmology=None,unit=dimensionless_unscaled):

		#Sanity check
		assert (cosmology is not None) or (comoving_distance is not None),"cosmology and comoving_distance cannot be both none!!"
		assert data.shape[0]==3,"A symmetric 2x2 matrix has 3 independent components!!"

		super(self.__class__,self).__init__(data,angle)
		self.redshift = redshift
		self.unit = unit

		#If a comoving distance is provided, use that; otherwise it needs to be computed from the astropy cosmology instance
		if comoving_distance is not None:		
			
			assert comoving_distance.unit.physical_type=="length"
			self.comoving_distance = comoving_distance
		
		else:
			self.comoving_distance = cosmology.comoving_distance(redshift)


#######################################################
###############RayTracer class#########################
#######################################################

class RayTracer(object):

	"""
	Class handler of ray tracing operations: it mainly computes the path corrections of light rays that travel through a set of gravitational lenses

	"""

	def __init__(self,lens_mesh_size=None):

		self.Nlenses = 0
		self.lens = list()
		self.distance = list()
		self.redshift = list()

		#If we know the size of the lens planes already we can compute, once and for all, the FFT meshgrid
		if lens_mesh_size is not None:
			self.lmesh = np.array(np.meshgrid(rfftfreq(lens_mesh_size),fftfreq(lens_mesh_size)))
		else:
			self.lmesh = None


	def addLens(self,lens_specification):

		"""
		Adds a gravitational lens to the ray tracer, either by putting in a lens plane, or by specifying the name of a file which contains the lens specifications

		:param lens_specification: specifications of the lens to add, either in tuple(filename,distance,redshift) or as a PotentialPlane instance
		:type lens specification: tuple or PotentialPlane instance

		"""

		#Sanity check
		assert type(lens_specification) in [tuple,PotentialPlane]

		#If specification is in tuple form, parse it
		if type(lens_specification)==tuple:

			filename,distance,redshift = lens_specification
			
			self.lens.append(filename)
			self.distance.append(distance)
			self.redshift.append(redshift)
			self.Nlenses += 1

		else:

			#Otherwise get the info from the PotentialPlane class
			self.lens.append(lens_specification)
			self.distance.append(lens_specification.comoving_distance)
			self.redshift.append(lens_specification.redshift)
			self.Nlenses += 1

		#If completed correctly, log info to the user
		logging.debug("Added lens at redshift {0:.3f}(comoving distance {1:.3f})".format(self.redshift[-1],self.distance[-1]))

	def randomRoll(self,seed=None):

		"""
		Randomly rolls all the lenses in the system along both axes

		:param seed: random seed with which to initialize the generator
		:type seed: int.

		"""

		if seed is not None:
			np.random.seed(seed)

		for lens in self.lens:
			lens.randomRoll(seed=None,lmesh=self.lmesh)


	def reorderLenses(self):

		"""
		Reorders the lenses in the ray tracer according to their comoving distance from the observer

		"""

		self.lens = [ lens for (redshift,lens) in sorted(zip(self.redshift,self.lens)) ]
		self.redshift.sort()
		self.distance.sort()


	##################################################################################################################################
	#####################This method solves the nonlinear lensing ODE#################################################################
	#############################(backward ray tracing)###############################################################################
	##################################################################################################################################


	def shoot(self,initial_positions,z=2.0,initial_deflection=None,precision="first",kind="positions",save_intermediate=False,compute_all_deflections=False,callback=None,**kwargs):

		"""
		Shots a bucket of light rays from the observer to the sources at redshift z (backward ray tracing), through the system of gravitational lenses, and computes the deflection statistics

		:param initial_positions: initial angular positions of the light ray bucket, according to the observer; if unitless, the positions are assumed to be in radians. initial_positions[0] is x, initial_positions[1] is y
		:type initial_positions: numpy array or quantity

		:param z: redshift of the sources; if an array is passed, a redshift must be specified for each ray, i.e. z.shape==initial_positions.shape[1:]
		:type z: float. or array

		:param initial_deflection: if not None, this is the initial deflection light rays undergo with respect to the line of sight (equivalent to specifying the first derivative IC on the lensing ODE); must have the same shape as initial_positions
		:type initial_deflection: numpy array or quantity

		:param precision: precision at which to compute weak lensing quantities, must be "first" for first order in the lensing potential, or "second" for added precision (not functional yet)
		:type precision: str.

		:param kind: what deflection statistics to compute; "positions" will calculate the ray deflections after they crossed the last lens, "jacobian" will compute the lensing jacobian matrix after the last lens, "shear" and "convergence" will compute the omonimous weak lensing statistics  
		:type kind: str.

		:param save_intermediate: save the intermediate positions of the rays too
		:type save_intermediate: bool.

		:param compute_all_deflections: if True, computes the gradients of the lensing potential at every pixel on the lens(might be overkill if Nrays<<Npixels); must be True if the computation is done with FFTs
		:type compute_all_deflections: bool.

		:param callback: if not None, this callback function is called on the current ray positions array at each step in the ray tracing; the current raytracing instance and the step number are passed as additional arguments, hence callback must match this signature
		:type callback: callable

		:param kwargs: the keyword arguments are passed to the callback if not None
		:type kwargs: dict.

		:returns: angular positions (or jacobians) of the light rays after the last lens crossing

		"""

		#Sanity check
		assert initial_positions.ndim>=2 and initial_positions.shape[0]==2,"initial positions shape must be (2,...)!"
		assert type(initial_positions)==quantity.Quantity and initial_positions.unit.physical_type=="angle"
		assert kind in ["positions","jacobians","shear","convergence"],"kind must be one in [positions,jacobians,shear,convergence]!"

		#Allocate arrays for the intermediate light ray positions and deflections

		if initial_deflection is None:
			current_positions = initial_positions.copy()
			current_deflection = np.zeros(initial_positions.shape) * initial_positions.unit
		else:
			assert initial_deflection.shape==initial_positions.shape
			current_deflection = initial_deflection.copy()
			current_positions = initial_positions + initial_deflection

		#If we want to trace jacobians, allocate also space for the jacobians
		if kind in ["jacobians","shear","convergence"]:

			#Initial condition for the jacobian is the identity
			current_jacobian = np.outer(np.array([1.0,0.0,0.0,1.0]),np.ones(initial_positions.shape[1:])).reshape((4,)+initial_positions.shape[1:])
			current_jacobian_deflection = np.zeros(current_jacobian.shape)

			#Useful to compute the product of the jacobian (2x2 matrix) with the shear matrix (2x2 symmetric matrix)
			dotter = np.zeros((4,3,4))
			dotter[(0,0,1,1,2,2,3,3),(0,2,0,2,2,1,2,1),(0,2,1,3,0,2,1,3)] = 1

		#Decide which is the last lens the light rays should cross
		if type(z)==np.ndarray:
			
			#Check that shapes correspond
			assert z.shape==initial_positions.shape[1:]

			#Check that redshift is not too high given the current lenses
			assert z.max()<self.redshift[-1],"Given the current lenses you can trace up to redshift {0:.2f}!".format(self.redshift[-1])

			#Compute the number of lenses that each ray should cross
			last_lens_ray = (z[None] > np.array(self.redshift).reshape((len(self.redshift),)+(1,)*len(z.shape))).argmin(0) - 1
			last_lens = last_lens_ray.max()
		
		else:
			
			#Check that redshift is not too high given the current lenses
			assert z<self.redshift[-1],"Given the current lenses you can trace up to redshift {0:.2f}!".format(self.redshift[-1])
			last_lens = (z>np.array(self.redshift)).argmin() - 1
		
		if kind=="positions" and save_intermediate:
			all_positions = np.zeros((last_lens+1,) + initial_positions.shape) * initial_positions.unit

		#The light rays positions at the k+1 th step are computed according to Xk+1 = Xk + Dk, where Dk is the deflection
		#To stabilize the solution numerically we compute the deflections as Dk+1 = (Ak-1)Dk + Ck*pk where pk is the deflection due to the potential gradient

		#Ordered references to the lenses
		distance = np.array([ d.to(Mpc).value for d in [0.0*Mpc] + self.distance ])
		redshift = np.array([0.0] + self.redshift)
		lens = self.lens

		#This is the main loop that goes through all the lenses
		for k in range(last_lens+1):

			#Load in the lens
			if type(lens[k])==PotentialPlane:
				current_lens = lens[k]
			elif type(lens[k])==str:
				current_lens = PotentialPlane.load(lens[k])
			else:
				raise TypeError("Lens format not recognized!")

			#Log
			logging.debug("Crossing lens {0} at redshift z={1:.2f}".format(k,current_lens.redshift))
			start = time.time()
			last_timestamp = start

			#Compute the deflection angles and log timestamp
			if compute_all_deflections:
				deflections = current_lens.deflectionAngles(lmesh=self.lmesh).getValues(current_positions[0],current_positions[1])
			else:
				deflections = current_lens.deflectionAngles(current_positions[0],current_positions[1])

			now = time.time()
			logging.debug("Retrieval of deflection angles from potential planes completed in {0:.3f}s".format(now-last_timestamp))
			last_timestamp = now

			#If we are tracing jacobians we need to retrieve the shear matrices too
			if kind in ["jacobians","convergence","shear"]:

				if compute_all_deflections:
					shear_tensors = current_lens.shearMatrix(lmesh=self.lmesh).getValues(current_positions[0],current_positions[1])
				else:
					shear_tensors = current_lens.shearMatrix(current_positions[0],current_positions[1])

				now = time.time()
				logging.debug("Shear matrices retrieved in {0:.3f}s".format(now-last_timestamp))
				last_timestamp = now
			
			#####################################################################################

			#Compute geometrical weight factors
			Ak = (distance[k+1] / distance[k+2]) * (1.0 + (distance[k+2] - distance[k+1])/(distance[k+1] - distance[k]))
			Ck = -1.0 * (distance[k+2] - distance[k+1]) / distance[k+2]

			#Compute the position on the next lens and log timestamp
			current_deflection *= (Ak-1) 
			now = time.time()
			logging.debug("Geometrical weight factors calculations and deflection scaling completed in {0:.3f}s".format(now-last_timestamp))
			last_timestamp = now

			#Add deflections and log timestamp
			current_deflection += Ck * deflections 
			now = time.time()
			logging.debug("Deflection angles computed in {0:.3f}s".format(now-last_timestamp))
			last_timestamp = now

			#If we are tracing jacobians we need to compute the matrix product with the shear matrix
			if kind in ["jacobians","convergence","shear"]:

				current_jacobian_deflection *= (Ak-1)

				#This is the part in which the products with the shear matrix are computed
				current_jacobian_deflection += Ck * (np.tensordot(dotter,current_jacobian,axes=([2],[0])) * shear_tensors).sum(1)
				
				now = time.time()
				logging.debug("Shear matrix products computed in {0:.3f}s".format(now-last_timestamp))
				last_timestamp = now

			###########################################################################################

			if type(z)==np.ndarray:
				
				current_positions[:,k<last_lens_ray] += current_deflection[:,k<last_lens_ray]
				current_positions[:,k==last_lens_ray] += current_deflection[:,k==last_lens_ray] * (z[None,k==last_lens_ray] - redshift[k+1]) / (redshift[k+2] - redshift[k+1])

				#We need to add the distortions to the jacobians too
				if kind in ["jacobians","convergence","shear"]:
					current_jacobian[:,k<last_lens_ray] += current_jacobian_deflection[:,k<last_lens_ray]
					current_jacobian[:,k==last_lens_ray] += current_jacobian_deflection[:,k==last_lens_ray] * (z[None,k==last_lens_ray] - redshift[k+1]) / (redshift[k+2] - redshift[k+1])

			else:
				
				if k<last_lens:
					current_positions += current_deflection
				else:
					current_positions += current_deflection * (z - redshift[k+1]) / (redshift[k+2] - redshift[k+1])

				#We need to add the distortions to the jacobians too
				if kind in ["jacobians","convergence","shear"]:

					if k<last_lens:
						current_jacobian += current_jacobian_deflection
					else:
						current_jacobian += current_jacobian_deflection * (z - redshift[k+1]) / (redshift[k+2] - redshift[k+1])

			now = time.time()
			logging.debug("Addition of deflections completed in {0:.3f}s".format(now-last_timestamp))
			last_timestamp = now

			#Save the intermediate positions if option was specified
			if kind=="positions" and save_intermediate:
				all_positions[k] = current_positions.copy()

			#Optionally, call the callback function on the current positions
			if callback is not None:
				callback(current_positions,self,k,**kwargs)

			#Log timestamp to cross lens
			now = time.time()
			logging.debug("Lens {0} crossed in {1:.3f}s".format(k,now-start))

		#Return the final positions of the light rays (or jacobians)
		if kind=="positions":
			
			if save_intermediate:
				return all_positions
			else:
				return current_positions

		else:

			#Different return types according to option (can compute convergence and shear directly)

			if kind=="convergence":
				return 1.0 - 0.5*(current_jacobian[0]+current_jacobian[3]) 
			
			elif kind=="shear":
				return np.array([0.5*(current_jacobian[3] - current_jacobian[0]),-0.5*(current_jacobian[1]+current_jacobian[2])])

			else:
				return current_jacobian

	#########################################################
	###########Direct calculation of the convergence#########
	#########################################################

	def convergenceDirect(self,initial_positions,z=2.0,save_intermediate=False):

		"""
		Computes the convergence directly integrating the lensing density along the unperturbed line of sight

		:param initial_positions: initial angular positions of the light ray bucket, according to the observer; if unitless, the positions are assumed to be in radians. initial_positions[0] is x, initial_positions[1] is y
		:type initial_positions: numpy array or quantity

		:param z: redshift of the sources
		:type z: float.

		:param save_intermediate: save the intermediate values of the convergence as successive lenses are crossed
		:type save_intermediate: bool.

		:returns: convergence values at each of the initial positions

		"""

		#Sanity check
		assert initial_positions.ndim>=2 and initial_positions.shape[0]==2,"initial positions shape must be (2,...)!"
		assert type(initial_positions)==quantity.Quantity and initial_positions.unit.physical_type=="angle"

		#Check that redshift is not too high given the current lenses
		assert z<self.redshift[-1],"Given the current lenses you can trace up to redshift {0:.2f}!".format(self.redshift[-1])
		last_lens = (z>np.array(self.redshift)).argmin() - 1

		if save_intermediate:
			all_convergence = np.zeros((last_lens+1,) + initial_positions.shape[1:])

		#Ordered references to the lenses
		distance = np.array([ d.to(Mpc).value for d in [0.0*Mpc] + self.distance ])
		redshift = np.array([0.0] + self.redshift)
		lens = self.lens

		#Timestamp
		now = time.time()
		last_timestamp = now

		#Loop that goes through the lenses
		current_convergence = np.zeros(initial_positions.shape[1:])
		for k in range(last_lens+1):

			#Start time for this lens
			start = time.time()

			#Load in the lens
			if type(lens[k])==PotentialPlane:
				current_lens = lens[k]
			elif type(lens[k])==str:
				current_lens = PotentialPlane.load(lens[k])
			else:
				raise TypeError("Lens format not recognized!")

			#Extract the density at the ray positions
			now = time.time()
			logging.debug("Extracting density values from lens {0} at redshift {1:2f}".format(k,current_lens.redshift))
			last_timestamp = now

			#Compute full density plane
			density = current_lens.density(initial_positions[0],initial_positions[1])

			#Timestamp
			now = time.time()
			logging.debug("Density values extracted in {0:.3f}s".format(now-last_timestamp))
			last_timestamp = now

			#Cumulate on the convergence
			if k<last_lens:
				current_convergence += 0.5 * density
			else:
				current_convergence += 0.5 * density * (z - redshift[k+1]) / (redshift[k+2] - redshift[k+1])

			now = time.time()
			logging.debug("Lens {0} crossed in {1:.3f}s".format(k,now-start))
			last_timestamp = now

			#Save the intermediate convergence values if option is enabled
			if save_intermediate:
				all_convergence[k] = current_convergence.copy()


		#Return to the user
		if save_intermediate:
			return all_convergence
		else:
			return current_convergence


	#########################################################
	############Forward ray tracing##########################
	#########################################################

	def shootForward(self,source_positions,z=2.0,save_intermediate=False,grid_resolution=512,interpolation="nearest"):

		"""
		Shoots a bucket of light rays from the source at redshift z to the observer at redshift 0 (forward ray tracing) and computes the according deflections using backward ray tracing plus a suitable interpolation scheme (KD Tree based)

		:param source_positions: angular positions of the unlensed sources
		:type source_positions: numpy array or quantity

		:param z: redshift of the sources
		:type z: float.

		:param save_intermediate: if True computes and saves the apparent image distortions after each lens is crossed (can be computationally expensive) 
		:type save_intermediate: bool.

		:param grid_resolution: the number of points on a side of the interpolation grid (must be choosen big enough according to the number of sources to resolve)
		:type grid_resolution: int. 

		:param interpolation: only "nearest" implemented so far
		:type interpolation: str.

		:returns: apparent positions of the sources as seen from the observer

		"""

		#First allocate the regular grid to use (must be fine enough to resolve the single ray distortions, since we use nearest neighbors interpolation for now)
		corner = source_positions.max(axis=tuple(range(1,len(source_positions.shape))))
		initial_grid = np.array(np.meshgrid(np.linspace(0.0,corner[0].value,grid_resolution),np.linspace(0.0,corner[1].value,grid_resolution))) * corner.unit
		initial_grid = initial_grid.reshape((2,)+(reduce(mul,initial_grid.shape[1:]),))

		now = time.time()
		last_timestamp = now
		
		#Perform the backwards ray tracing
		final_grid = self.shoot(initial_grid,z=z,save_intermediate=save_intermediate)

		now = time.time()
		logging.debug("Ray tracing in {0:.3f}s".format(now-last_timestamp))
		last_timestamp = now

		if save_intermediate:

			#If this option is enabled the full evolution of the distortions (after each lens is crossed) is computed

			#Allocate space for the apparent positions
			apparent_positions = np.zeros((final_grid.shape[0],) + source_positions.shape) * corner.unit

			for n in range(final_grid.shape[0]):

				#Next build the KD tree for nearest neighbors interpolation
				tree = KDTree(final_grid[n].transpose())

				now = time.time()
				logging.debug("KDTree built in {0:.3f}s".format(now-last_timestamp))
				last_timestamp = now

				#Query the tree and retrieve the apparent positions
				distances,apparent_position_index = tree.query(source_positions.reshape((2,)+(reduce(mul,source_positions.shape[1:]),)).transpose())

				now = time.time()
				logging.debug("Tree query completed in {0:.3f}s".format(now-last_timestamp))
				last_timestamp = now

				apparent_positions[n] = initial_grid[:,apparent_position_index].reshape(source_positions.shape).copy()


		else:

			#Next build the KD tree for nearest neighbors interpolation
			tree = KDTree(final_grid.transpose())

			now = time.time()
			logging.debug("KDTree built in {0:.3f}s".format(now-last_timestamp))
			last_timestamp = now

			#Query the tree and retrieve the apparent positions
			distances,apparent_position_index = tree.query(source_positions.reshape((2,)+(reduce(mul,source_positions.shape[1:]),)).transpose())

			now = time.time()
			logging.debug("Tree query completed in {0:.3f}s".format(now-last_timestamp))
			last_timestamp = now

			apparent_positions = initial_grid[:,apparent_position_index].reshape(source_positions.shape)

		#Return the measured apparent distances
		return apparent_positions



	#########################################################
	#######################Graphics##########################
	#########################################################

	def reset(self):

		"""
		Resets the RayTracer plotting engine
		"""

		self._x_pos = list()
		self._y_pos = list()

	def displayRays(self,initial_positions,z=2.0,projection="2d",fig=None,ax=None):

		"""
		Uses matplotlib to display a visual of the lens system and of the deflection that the light rays which traverse it undergo

		param initial_positions: initial angular positions of the light ray bucket, according to the observer; if unitless, the positions are assumed to be in radians. initial_positions[0] is x, initial_positions[1] is y
		:type initial_positions: numpy array or quantity

		:param z: redshift of the sources; if an array is passes, a redshift must be specified for each ray, i.e. z.shape==initial_positions.shape[1:]
		:type z: float. or array

		:param projection: can be "2d" for the projections of the ray positions along x and y or "3d" for the full 3d view
		:type projection: str.

		:param fig: figure object that owns the plot
		:type fig: matplotlib figure

		"""

		if matplotlib is None:
			raise ValueError("matplotlib not found!")

		#Instantiate axes
		if (fig is None) or (ax is None):

			if projection=="2d":
				self.fig,self.ax = plt.subplots(1,2,figsize=(16,8),subplot_kw={"axisbg":"black"})
			elif projection=="3d":
				pass
			else:
				raise ValueError("projection must be either 2d or 3d!")

		else:
			self.fig = fig
			self.ax = ax

		#Compute the positions of the light rays across the lens system
		initial_positions = initial_positions.reshape((2,)+(reduce(mul,initial_positions.shape[1:]),))
		pos = self.shoot(initial_positions,z,save_intermediate=True)

		#Add a 0 position corresponding to the observer, and the initial positions too
		pos = np.concatenate((np.zeros((1,) + pos.shape[1:]),initial_positions[None,:].value,pos.value)) * pos.unit

		#Construct an array with the comoving distance of the lenses
		distance_lenses = np.array([ d.to(Mpc).value for d in [0.0 * Mpc] + self.distance[:pos.shape[0]] ]) * Mpc
		distance = distance_lenses.copy()
		redshift = np.array([0.0]+self.redshift[:pos.shape[0]])

		#Correct the last position if the last redshift falls between two lenses
		distance[-1] = distance[-2] + (distance[-1] - distance[-2]) * (z - redshift[-2]) / (redshift[-1] - redshift[-2])

		if projection=="2d":

			if not (hasattr(self,"_x_pos") and hasattr(self,"_y_pos")):
				self.reset()

			#Plot the x,y positions
			if len(self._x_pos)==0 and len(self._y_pos)==0:
			
				for nray in range(pos.shape[2]):

					self._x_pos.append(self.ax[0].plot(distance[:pos.shape[0]],distance[:pos.shape[0]]*pos[:,0,nray].to(rad),color="yellow")[0])
					self._y_pos.append(self.ax[1].plot(distance[:pos.shape[0]],distance[:pos.shape[0]]*pos[:,1,nray].to(rad),color="yellow")[0])
				
				self.ax[0].set_xlabel(r"$w$({0})".format(distance.unit.to_string()))
				self.ax[0].set_ylabel(r"$x$({0})".format(distance.unit.to_string()))
				self.ax[1].set_xlabel(r"$w$({0})".format(distance.unit.to_string()))
				self.ax[1].set_ylabel(r"$y$({0})".format(distance.unit.to_string()))

				#Plot the lenses too
				for d in distance_lenses:
					for i in range(2):
						min = distance[-1]*pos.to(rad).value.min()
						max = distance[-1]*pos.to(rad).value.max()
						self.ax[i].plot(d*np.ones(100),np.linspace(min,max,100),color="red")


			else:

				for nray in range(pos.shape[2]):
					self._x_pos[nray].set_xdata(distance[:pos.shape[0]])
					self._y_pos[nray].set_xdata(distance[:pos.shape[0]])
					self._x_pos[nray].set_ydata(distance[:pos.shape[0]]*pos[:,0,nray].to(rad))
					self._y_pos[nray].set_ydata(distance[:pos.shape[0]]*pos[:,1,nray].to(rad))


				#Plot the lenses too
				for d in distance_lenses:
					for i in range(2):
						min = distance[-1]*pos.to(rad).value.min()
						max = distance[-1]*pos.to(rad).value.max()
						self.ax[i].plot(d*np.ones(100),np.linspace(min,max,100),color="red")


			#Adjust the x,y axis range
			self.ax[0].set_xlim(0,self.distance[-1].value)
			self.ax[1].set_xlim(0,self.distance[-1].value)

		else:
			pass






		

