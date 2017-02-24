from __future__ import division

import re

import numpy as np
import astropy.units as u
from astropy.cosmology import w0waCDM

#Try to import the FITSIO library for optimal FITS images reading
try:
	from fitsio import FITS as fitsio
	fitsio = fitsio
	from astropy.io import fits

except ImportError:
	from astropy.io import fits
	fitsio = None

####################################################################
#######################FITS format##################################
####################################################################

#Header
def readFITSHeader(filename):
	with fits.open(filename) as fp:
		return fp[0].header

#Read
def readFITS(cls,filename,init_cosmology=True):

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
		hubble = header["H0"] * (u.km/(u.s*u.Mpc))
		h = header["h"]
	except:
		hubble = header["H_0"] * (u.km/(u.s*u.Mpc))
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
	comoving_distance = (header["CHI"] / h) * u.Mpc

	if "SIDE" in header.keys():
		angle = header["SIDE"] * u.Mpc / h
	elif "ANGLE" in header.keys():
		angle = header["ANGLE"] * u.deg
	else:
		angle = ((header["RES_X"] * header["NAXIS1"] / header["CHI"]) * u.rad).to(u.deg)

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
		unit = getattr(u,name)
		if exponent is not None:
			unit **= exponent
	except AttributeError:
		unit = u.dimensionless_unscaled
	except (ValueError,KeyError):
		unit = u.rad**2

	#Instantiate the new PotentialPlane instance
	if fitsio is not None:

		if len(hdu)==1:
			new_plane = cls(hdu[0].read(),angle=angle,redshift=redshift,comoving_distance=comoving_distance,cosmology=cosmology,unit=unit,num_particles=num_particles,filename=filename)
		else:
			new_plane = cls(hdu[1].read() + 1.0j*hdu[1].read(),angle=angle,redshift=redshift,comoving_distance=comoving_distance,cosmology=cosmology,unit=unit,num_particles=num_particles,filename=filename)

	else:
			
		if len(hdu)==1:
			new_plane = cls(hdu[0].data.astype(np.float64),angle=angle,redshift=redshift,comoving_distance=comoving_distance,cosmology=cosmology,unit=unit,num_particles=num_particles,filename=filename)
		else:
			new_plane = cls((hdu[0].data + 1.0j*hdu[1].data).astype(np.complex128),angle=angle,redshift=redshift,comoving_distance=comoving_distance,cosmology=cosmology,unit=unit,num_particles=num_particles,filename=filename)

	#Close the FITS file and return
	hdu.close()
	return new_plane



#Write
def saveFITS(self,filename,double_precision):

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
		raise ValueError("Space must either be real of Fourier!")


	#Generate a header
	hdu.header["H0"] = (self.cosmology.H0.to(u.km/(u.s*u.Mpc)).value,"Hubble constant in km/s*Mpc")
	hdu.header["h"] = (self.cosmology.h,"Dimensionless Hubble constant")
	hdu.header["OMEGA_M"] = (self.cosmology.Om0,"Dark Matter density")
	hdu.header["OMEGA_L"] = (self.cosmology.Ode0,"Dark Energy density")
	hdu.header["W0"] = (self.cosmology.w0,"Dark Energy equation of state")
	hdu.header["WA"] = (self.cosmology.wa,"Dark Energy running equation of state")

	hdu.header["Z"] = (self.redshift,"Redshift of the lens plane")
	hdu.header["CHI"] = (hdu.header["h"] * self.comoving_distance.to(u.Mpc).value,"Comoving distance in Mpc/h")

	if self.side_angle.unit.physical_type=="angle":
		hdu.header["ANGLE"] = (self.side_angle.to(u.deg).value,"Side angle in degrees")
	elif self.side_angle.unit.physical_type=="length":
		hdu.header["SIDE"] = (self.side_angle.to(u.Mpc).value*self.cosmology.h,"Side length in Mpc/h")

	hdu.header["NPART"] = (float(self.num_particles),"Number of particles on the plane")
	hdu.header["UNIT"] = (self.unit.to_string(),"Pixel value unit") 

	#Save the plane
	if self.space=="real":
		hdulist = fits.HDUList([hdu])
	else:
		hdulist = fits.HDUList([hdu,hdu1])

	hdulist.writeto(filename,overwrite=True)

########################################################################################################################################