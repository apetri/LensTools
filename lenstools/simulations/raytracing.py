from ..convergence import Spin0
from ..shear import Spin1

import numpy as np
from astropy.units import km,s,Mpc,rad,deg
from astropy.io import fits

###########################################################
###############PotentialPlane class########################
###########################################################

class PotentialPlane(Spin0):

	"""
	Class handler of a lens potential plane, inherits from the parent Spin0 class; additionally it defines redshift and comoving distance attributes that are needed for ray tracing operations

	"""

	def __init__(self,data,angle,redshift,cosmology,unit=rad**2,masked=False):

		super(self.__class__,self).__init__(data,angle,masked)
		self.redshift = redshift
		self.comoving_distance = cosmology.comoving_distance(redshift)
		self.cosmology = cosmology
		self.unit = unit


	def save(self,filename,format="fits"):

		"""
		Saves the potential plane to an external file, of which the format can be specified (only fits implemented so far)

		"""

		#Create the hdu
		hdu = fits.PrimaryHDU(self.data)

		#Generate a header
		hdu.header["H0"] = self.cosmology.H0.to(km/(s*Mpc)).value
		hdu.header["h"] = self.cosmology.h
		hdu.header["OMEGA_M"] = self.cosmology.Om0
		hdu.header["OMEGA_L"] = self.cosmology.Ode0
		hdu.header["W0"] = self.cosmology.w0
		hdu.header["WA"] = self.cosmology.wa

		hdu.header["Z"] = (self.redshift,"Redshift of the lens plane")
		hdu.header["CHI"] = (hdu.header["h"] * self.comoving_distance.to(Mpc).value,"Comoving distance in Mpc/h")

		hdu.header["ANGLE"] = (self.angle.to(deg).value,"Side angle in degrees")

		#Save the plane
		hdulist = fits.HDUList([hdu])
		hdulist.writeto(filename)


	
	def deflectionAngles():

		"""
		Computes the deflection angles for the given lensing potential by taking the gradient of the potential map

		"""

		#Compute the gradient of the potential map
		deflection_x,deflection_y = self.gradient() / self.resolution.to(self.unit**0.5).value

		#Return the DeflectionPlane instance
		return DeflectionPlane(np.array(deflection_x,deflection_y),angle=self.angle,redshift=self.redshift,cosmology=self.cosmology,unit=self.unit**0.5)




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
	Docstring
	"""