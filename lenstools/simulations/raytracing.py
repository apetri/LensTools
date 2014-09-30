from ..convergence import Spin0
from ..shear import Spin1

import logging

import numpy as np

#FFT engines
from numpy.fft import fftfreq,rfft2,irfft2

try:
	from numpy.fft import rfftfreq
except ImportError:
	from ..utils import rfftfreq

from astropy.cosmology import w0waCDM
from astropy.units import km,s,Mpc,rad,deg
from astropy.io import fits

###########################################################
###############PotentialPlane class########################
###########################################################

class PotentialPlane(Spin0):

	"""
	Class handler of a lens potential plane, inherits from the parent Spin0 class; additionally it defines redshift and comoving distance attributes that are needed for ray tracing operations

	"""

	def __init__(self,data,angle,redshift,cosmology,unit=rad**2,space="real",masked=False):

		super(self.__class__,self).__init__(data,angle,masked)
		self.redshift = redshift
		self.comoving_distance = cosmology.comoving_distance(redshift)
		self.cosmology = cosmology
		self.unit = unit
		self.space = space


	def save(self,filename,format="fits"):

		"""
		Saves the potential plane to an external file, of which the format can be specified (only fits implemented so far)

		:param filename: name of the file on which to save the plane
		:type filename: str.

		:param format: format of the file, only FITS implemented so far
		:type format: str.

		"""

		if format=="fits":
		
			#Create the hdu
			hdu = fits.PrimaryHDU(self.data)

			#Generate a header
			hdu.header["H0"] = (self.cosmology.H0.to(km/(s*Mpc)).value,"Hubble constant in km/s*Mpc")
			hdu.header["h"] = (self.cosmology.h,"Dimensionless Hubble constant")
			hdu.header["OMEGA_M"] = (self.cosmology.Om0,"Dark Matter density")
			hdu.header["OMEGA_L"] = (self.cosmology.Ode0,"Dark Energy density")
			hdu.header["W0"] = (self.cosmology.w0,"Dark Energy equation of state")
			hdu.header["WA"] = (self.cosmology.wa,"Dark Energy running equation of state")

			hdu.header["Z"] = (self.redshift,"Redshift of the lens plane")
			hdu.header["CHI"] = (hdu.header["h"] * self.comoving_distance.to(Mpc).value,"Comoving distance in Mpc/h")

			hdu.header["ANGLE"] = (self.side_angle.to(deg).value,"Side angle in degrees")

			#Save the plane
			hdulist = fits.HDUList([hdu])
			hdulist.writeto(filename,clobber=True)

		else:
			raise ValueError("Format {0} not implemented yet!!".format(format))


	@classmethod
	def load(cls,filename,format="fits"):

		"""
		Loads the potential plane from an external file, of which the format can be specified (only fits implemented so far)

		:param filename: name of the file from which to load the plane
		:type filename: str.

		:param format: format of the file, only FITS implemented so far
		:type format: str.

		:returns: PotentialPlane instance that wraps the data contained in the file

		"""

		if format=="fits":

			#Read the FITS file with the plane information
			hdu = fits.open(filename)

			#Retrieve the info from the header
			hubble = hdu[0].header["H0"] * (km/(s*Mpc))
			Om0 = hdu[0].header["OMEGA_M"]
			Ode0 = hdu[0].header["OMEGA_L"]
			w0 = hdu[0].header["W0"]
			wa = hdu[0].header["WA"]
			redshift = hdu[0].header["Z"]
			angle = hdu[0].header["ANGLE"] * deg

			#Build the cosmology object
			cosmology = w0waCDM(H0=hubble,Om0=Om0,Ode0=Ode0,w0=w0,wa=wa)

			#Instantiate the new PotentialPlane instance
			return cls(hdu[0].data.astype(np.float64),angle=angle,redshift=redshift,cosmology=cosmology,unit=rad**2)

		else:
			raise ValueError("Format {0} not implemented yet!!".format(format))


	
	def deflectionAngles(self):

		"""
		Computes the deflection angles for the given lensing potential by taking the gradient of the potential map; it is also possible to proceed with FFTs

		"""

		#Compute the gradient of the potential map

		if self.space=="real":
			
			deflection_x,deflection_y = self.gradient()
		
		elif self.space=="fourier":

			#Compute deflections in fourier space
			ly,lx = np.meshgrid(fftfreq(self.data.shape[0]),rfftfreq(self.data.shape[0]),indexing="ij")
			ft_deflection_x = 1.0j * self.data * lx 
			ft_deflection_y = 1.0j * self.data * ly

			#Go back in real space
			deflection_x = irfft2(ft_deflection_x)
			deflection_y = irfft2(ft_deflection_y)

		else:
			raise ValueError("space must be either real or fourier!")

		#Return the DeflectionPlane instance
		return DeflectionPlane(np.array([deflection_x,deflection_y])/self.resolution.to(self.unit**0.5).value,angle=self.side_angle,redshift=self.redshift,cosmology=self.cosmology,unit=self.unit**0.5)


	def toReal(self):

		"""
		Switches to real space

		"""

		assert self.space=="fourier"
		self.data = irfft2(self.data)
		self.space = "real"

	
	def toFourier(self):

		"""
		Switches to Fourier space

		"""

		assert self.space=="real"
		self.data = rfft2(self.data)
		self.space="fourier"


	def density(self):

		"""
		Computes the projected density fluctuation by taking the laplacian of the potential; useful to check if the potential is reasonable

		:returns: Spin0 instance with the density fluctuation data 

		"""

		assert self.space=="real","You really want to do this operation in fourier space?"

		#Compute the laplacian
		hessian_xx,hessian_yy,hessian_xy = self.hessian()

		#The density is twice the trace of the hessian
		return Spin0(2.0*(hessian_xx + hessian_yy)/(self.resolution**2).to(self.unit).value,angle=self.side_angle)




#############################################################
################DeflectionPlane class########################
#############################################################

class DeflectionPlane(Spin1):

	"""
	Class handler of a lens deflection plane, inherits from the parent Spin1 class and holds the values of the deflection angles of the light rays that cross a potential plane

	"""

	def __init__(self,data,angle,redshift,cosmology,unit=rad):

		super(self.__class__,self).__init__(data,angle)
		self.redshift = redshift
		self.comoving_distance = cosmology.comoving_distance(redshift)
		self.cosmology = cosmology
		self.unit = unit



#######################################################
###############RayTracer class#########################
#######################################################

class RayTracer(object):

	"""
	Class handler of ray tracing operations: it mainly computes the path corrections of light rays that travel through a set of gravitational lenses

	"""

	def __init__(self):

		self.Nlenses = 0
		self.lens = list()
		self.distance = list()
		self.redshift = list()


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


	def shoot(self,initial_positions,z=2.0,precision="first",kind="positions"):

		"""
		Shots a bucket of light rays from the observer to the sources at redshift z, through the system of gravitational lenses, and computes the deflection statistics

		:param initial_positions: initial angular positions of the light ray bucket, according to the observer; if unitless, the positions are assumed to be in radians. initial_posiions[0] is x, initial_positions[1] is y
		:type initial_positions: numpy array or quantity

		:param z: redshift of the sources
		:type z: float.

		:param precision: precision at which to compute weak lensing quantities, must be "first" for first order in the lensing potential, or "second" for added precision
		:type precision: str.

		:param kind: what deflection statistics to compute; "positions" will calculate the ray deflections after they crossed the last lens, "jacobian" will compute the lensing jacobian matrix after the last lens, "shear" and "convergence" will compute the omonimous weak lensing statistics  
		:type kind: str.

		"""

		return None

