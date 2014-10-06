from ..convergence import Spin0,ConvergenceMap
from ..shear import Spin1,Spin2,ShearMap

import time
import logging

import numpy as np

#FFT engines
from numpy.fft import fftfreq,rfft2,irfft2

try:
	from numpy.fft import rfftfreq
except ImportError:
	from ..utils import rfftfreq

from astropy.cosmology import w0waCDM
from astropy.units import km,s,Mpc,rad,deg,dimensionless_unscaled,quantity
from astropy.io import fits

###########################################################
###############PotentialPlane class########################
###########################################################

class PotentialPlane(Spin0):

	"""
	Class handler of a lens potential plane, inherits from the parent Spin0 class; additionally it defines redshift and comoving distance attributes that are needed for ray tracing operations

	"""

	def __init__(self,data,angle,redshift,cosmology=None,comoving_distance=None,unit=rad**2,masked=False):

		#Sanity check
		assert (cosmology is not None) or (comoving_distance is not None),"cosmology and comoving_distance cannot be both None!!"

		super(self.__class__,self).__init__(data,angle,masked)
		self.redshift = redshift
		self.cosmology = cosmology
		self.unit = unit

		#If a comoving distance is provided, use that; otherwise it needs to be computed from the astropy cosmology instance
		if comoving_distance is not None:		
			
			assert comoving_distance.unit.physical_type=="length"
			self.comoving_distance = comoving_distance
		
		else:
			self.comoving_distance = cosmology.comoving_distance(redshift)

		if data.dtype==np.float:
			self.space = "real"
		elif data.dtype==np.complex:
			self.space = "fourier"
		else:
			raise TypeError("data type not supported!")


	def save(self,filename,format="fits"):

		"""
		Saves the potential plane to an external file, of which the format can be specified (only fits implemented so far)

		:param filename: name of the file on which to save the plane
		:type filename: str.

		:param format: format of the file, only FITS implemented so far
		:type format: str.

		"""

		if format=="fits":

			#A cosmology instance should be available in order to save in FITS format
			assert self.cosmology is not None
		
			#Create the hdu
			if self.space=="real":
				
				hdu = fits.PrimaryHDU(self.data)
			
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

			hdu.header["ANGLE"] = (self.side_angle.to(deg).value,"Side angle in degrees")

			#Save the plane
			if self.space=="real":
				hdulist = fits.HDUList([hdu])
			else:
				hdulist = fits.HDUList([hdu,hdu1])

			hdulist.writeto(filename,clobber=True)

		else:
			raise ValueError("Format {0} not implemented yet!!".format(format))


	@classmethod
	def load(cls,filename,format="fits",init_cosmology=True):

		"""
		Loads the potential plane from an external file, of which the format can be specified (only fits implemented so far)

		:param filename: name of the file from which to load the plane
		:type filename: str.

		:param format: format of the file, only FITS implemented so far
		:type format: str.

		:param init_cosmology: if True, instantiates the cosmology attribute of the PotentialPlane
		:type init_cosmology: bool.

		:returns: PotentialPlane instance that wraps the data contained in the file

		"""

		if format=="fits":

			#Read the FITS file with the plane information (if there are two HDU's the second one is the imaginary part)
			hdu = fits.open(filename)
			if len(hdu)>2:
				raise ValueError("There are more than 2 HDUs, file format unknown")

			header = hdu[0].header

			#Retrieve the info from the header
			hubble = header["H0"] * (km/(s*Mpc))
			h = header["h"]
			Om0 = header["OMEGA_M"]
			Ode0 = header["OMEGA_L"]
			w0 = header["W0"]
			wa = header["WA"]
			redshift = header["Z"]
			angle = header["ANGLE"] * deg
			comoving_distance = (header["CHI"] / h) * Mpc

			#Build the cosmology object if options directs
			if init_cosmology:
				cosmology = w0waCDM(H0=hubble,Om0=Om0,Ode0=Ode0,w0=w0,wa=wa)
			else:
				cosmology = None

			#Instantiate the new PotentialPlane instance
			if len(hdu)==1:
				return cls(hdu[0].data.astype(np.float64),angle=angle,redshift=redshift,comoving_distance=comoving_distance,cosmology=cosmology,unit=rad**2)
			else:
				return cls((hdu[0].data + 1.0j*hdu[1].data).astype(np.complex128),angle=angle,redshift=redshift,comoving_distance=comoving_distance,cosmology=cosmology,unit=rad**2)

			#Close the FITS file
			hdu.close()

		else:
			raise ValueError("Format {0} not implemented yet!!".format(format))


	def randomRoll(self,seed=None):

		"""
		Randomly shifts the plane along its axes, enforcing periodic boundary conditions

		:param seed: random seed with which to initialize the generator
		:type seed: int.

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
			ly,lx = np.meshgrid(fftfreq(self.data.shape[0]),rfftfreq(self.data.shape[0]),indexing="ij")

			#Timestamp
			now = time.time()
			logging.debug("l meshgrid initialized in {0:.3f}s".format(now-last_timestamp))
			last_timestamp = now 

			self.data *= np.exp(2.0j*np.pi*(lx*np.random.randint(0,self.data.shape[0]) + ly*np.random.randint(0,self.data.shape[0])))

			#Timestamp
			now = time.time()
			logging.debug("Phase multiplication completed in {0:.3f}s".format(now-last_timestamp))
			last_timestamp = now 


		else:
			raise ValueError("space must be either real or fourier!")

	
	def deflectionAngles(self):

		"""
		Computes the deflection angles for the given lensing potential by taking the gradient of the potential map; it is also possible to proceed with FFTs

		"""

		now = time.time()
		last_timestamp = now

		if self.space=="real":
			
			#Compute the gradient of the potential map
			deflection_x,deflection_y = self.gradient()
		
		elif self.space=="fourier":

			#Compute deflections in fourier space
			ly,lx = np.meshgrid(fftfreq(self.data.shape[0]),rfftfreq(self.data.shape[0]),indexing="ij")

			#Timestamp
			now = time.time()
			logging.debug("l meshgrid initialized in {0:.3f}s".format(now-last_timestamp))
			last_timestamp = now 

			ft_deflection_x = 2.0*np.pi*1.0j * self.data * lx 
			ft_deflection_y = 2.0*np.pi*1.0j * self.data * ly

			#Timestamp
			now = time.time()
			logging.debug("Phase multiplications completed in {0:.3f}s".format(now-last_timestamp))
			last_timestamp = now 

			#Go back in real space
			deflection_x = irfft2(ft_deflection_x)
			deflection_y = irfft2(ft_deflection_y)

			#Timestamp
			now = time.time()
			logging.debug("Inverse FFTs completed in {0:.3f}s".format(now-last_timestamp))
			last_timestamp = now 

		else:
			raise ValueError("space must be either real or fourier!")

		#Return the DeflectionPlane instance
		return DeflectionPlane(np.array([deflection_x,deflection_y])/self.resolution.to(self.unit**0.5).value,angle=self.side_angle,redshift=self.redshift,comoving_distance=self.comoving_distance,cosmology=self.cosmology,unit=self.unit**0.5)


	def shearMatrix(self):

		"""
		Computes the shear matrix for the given lensing potential; it is also possible to proceed with FFTs

		"""

		if self.space=="real":
			
			#Compute the second derivatives
			tensor_xx,tensor_yy,tensor_xy = self.hessian()

		elif self.space=="fourier":

			#Compute deflections in fourier space
			ly,lx = np.meshgrid(fftfreq(self.data.shape[0]),rfftfreq(self.data.shape[0]),indexing="ij")
			ft_tensor_xx = -(2.0*np.pi)**2 * self.data * (lx**2)
			ft_tensor_xy = -(2.0*np.pi)**2 * self.data * (lx*ly)
			ft_tensor_yy = -(2.0*np.pi)**2 * self.data * (ly**2)

			#Go back in real space
			tensor_xx = irfft2(ft_tensor_xx)
			tensor_xy = irfft2(ft_tensor_xy)
			tensor_yy = irfft2(ft_tensor_yy)

		else:
			raise ValueError("space must be either real or fourier!")


		#Return the ShearTensorPlane instance
		return ShearTensorPlane(np.array([tensor_xx,tensor_yy,tensor_xy])/self.resolution.to(self.unit**0.5).value**2,angle=self.side_angle,redshift=self.redshift,comoving_distance=self.comoving_distance,cosmology=self.cosmology,unit=dimensionless_unscaled)


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


	def density(self):

		"""
		Computes the projected density fluctuation by taking the laplacian of the potential; useful to check if the potential is reasonable

		:returns: Spin0 instance with the density fluctuation data 

		"""

		#Compute the laplacian
		if self.space=="real":			
			
			hessian_xx,hessian_yy,hessian_xy = self.hessian()
			laplacian = hessian_xx + hessian_yy

		elif self.space=="fourier":

			ly,lx = np.meshgrid(fftfreq(self.data.shape[0]),rfftfreq(self.data.shape[0]),indexing="ij")
			ft_laplacian = -1.0 * (2.0*np.pi)**2 * (lx**2 + ly**2) * self.data
			laplacian = irfft2(ft_laplacian) 			

		else:
			raise ValueError("space must be either real or fourier!")

		#The density is twice the trace of the hessian
		return Spin0(2.0*laplacian/(self.resolution**2).to(self.unit).value,angle=self.side_angle)




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

		return self.gradient() / self.resolution.to(self.unit).value


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

	def randomRoll(self,seed=None):

		"""
		Randomly rolls all the lenses in the system along both axes

		:param seed: random seed with which to initialize the generator
		:type seed: int.

		"""

		if seed is not None:
			np.random.seed(seed)

		for lens in self.lens:
			lens.randomRoll(seed=None)


	def reorderLenses(self):

		"""
		Reorders the lenses in the ray tracer according to their comoving distance from the observer

		"""

		self.lens = [ lens for (redshift,lens) in sorted(zip(self.redshift,self.lens)) ]
		self.redshift.sort()
		self.distance.sort()


	##################################################################################################################################
	#####################This method solves the nonlinear lensing ODE#################################################################
	##################################################################################################################################


	def shoot(self,initial_positions,z=2.0,precision="first",kind="positions",save_intermediate=False):

		"""
		Shots a bucket of light rays from the observer to the sources at redshift z, through the system of gravitational lenses, and computes the deflection statistics

		:param initial_positions: initial angular positions of the light ray bucket, according to the observer; if unitless, the positions are assumed to be in radians. initial_posiions[0] is x, initial_positions[1] is y
		:type initial_positions: numpy array or quantity

		:param z: redshift of the sources; if an array is passes, a redshift must be specified for each ray, i.e. z.shape==initial_positions.shape[1:]
		:type z: float. or array

		:param precision: precision at which to compute weak lensing quantities, must be "first" for first order in the lensing potential, or "second" for added precision
		:type precision: str.

		:param kind: what deflection statistics to compute; "positions" will calculate the ray deflections after they crossed the last lens, "jacobian" will compute the lensing jacobian matrix after the last lens, "shear" and "convergence" will compute the omonimous weak lensing statistics  
		:type kind: str.

		:param save_intermediate: save the intermediate positions of the rays too
		:type save_intermediate: bool.

		"""

		#Sanity check
		assert initial_positions.ndim>=2 and initial_positions.shape[0]==2,"initial positions shape must be (2,...)!"
		assert type(initial_positions)==quantity.Quantity and initial_positions.unit.physical_type=="angle"
		assert kind in ["positions","jacobians","shear","convergence"],"kind must be one in [positions,jacobians,shear,convergence]!"

		#Allocate arrays for the intermediate light ray positions
		current_positions = initial_positions.copy()
		current_deflection = np.zeros(initial_positions.shape) * initial_positions.unit

		#Decide which is the last lens the light rays should cross
		if type(z)==np.ndarray:
			
			#Check that shapes correspond
			assert z.shape==initial_positions.shape[1:]

			#Compute the number of lenses that each ray should cross
			last_lens_ray = (z[None] > np.array(self.redshift).reshape((len(self.redshift),)+(1,)*len(z.shape))).argmin(0) - 1
			last_lens = last_lens_ray.max()
		
		else:
			last_lens = (z>np.array(self.redshift)).argmin() - 1
		
		if save_intermediate:
			all_positions = np.zeros((last_lens,) + initial_positions.shape) * initial_positions.unit

		#Allocate a new array of dimensionless distances
		distance = np.array([ d.to(Mpc).value for d in [0.0*Mpc] + self.distance ])

		#The light rays positions at the k+1 th step are computed according to Xk+1 = Xk + Dk, where Dk is the deflection
		#To stabilize the solution numerically we compute the deflections as Dk+1 = (Ak-1)Dk + Ck*pk where pk is the deflection due to the potential gradient

		#This is the main loop that goes through all the lenses
		for k in range(last_lens):

			#Load in the lens
			if type(self.lens[k])==PotentialPlane:
				current_lens = self.lens[k]
			elif type(self.lens[k])==str:
				current_lens = PotentialPlane.load(self.lens[k])
			else:
				raise TypeError("Lens format not recognized!")

			#Log
			logging.debug("Crossing lens {0} at redshift z={1:.2f}".format(k,current_lens.redshift))
			start = time.time()
			last_timestamp = start

			#Compute the deflection angles and log timestamp
			deflections = current_lens.deflectionAngles()
			now = time.time()
			logging.debug("Deflection angles computed in {0:.3f}s".format(now-last_timestamp))
			last_timestamp = now

			#Compute geometrical weight factors
			Ak = (distance[k+1] / distance[k+2]) * (1.0 + (distance[k+2] - distance[k+1])/(distance[k+1] - distance[k]))
			Ck = -1.0 * (distance[k+2] - distance[k+1]) / distance[k+2]

			#Compute the position on the next lens and log timestamp
			current_deflection *= (Ak-1) 
			now = time.time()
			logging.debug("Geometrical weight factors calculations and deflection scaling completed in {0:.3f}s".format(now-last_timestamp))
			last_timestamp = now

			#Add deflections and log timestamp
			current_deflection += Ck * deflections.getValues(current_positions[0],current_positions[1]) * deflections.unit
			now = time.time()
			logging.debug("Retrieval of deflection angles from potential planes completed in {0:.3f}s".format(now-last_timestamp))
			last_timestamp = now


			if type(z)==np.ndarray:
				current_positions[:,k<last_lens_ray] += current_deflection[:,k<last_lens_ray]
			else:
				current_positions += current_deflection

			now = time.time()
			logging.debug("Addition of deflections completed in {0:.3f}s".format(now-last_timestamp))
			last_timestamp = now

			#Save the intermediate positions if option was specified
			if save_intermediate:
				all_positions[k] = current_positions.copy()

			#Log timestamp to cross lens
			now = time.time()
			logging.debug("Lens {0} crossed in {1:.3f}s".format(k,now-start))

		#Return the final positions of the light rays
		if save_intermediate:
			return all_positions
		else:
			return current_positions




		

