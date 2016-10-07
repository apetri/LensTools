import numpy as np
from astropy.io import fits
import astropy.units as u

#Load scalar FITS image
def loadFITS(cls,filename):

	with fits.open(filename) as hdu:
		data = hdu[0].data
		angle = hdu[0].header["ANGLE"] * u.deg
		return cls(data,angle)

#Write FITS scalar image
def saveFITS(self,filename,double_precision):
		
	#Create the hdu
	if double_precision:
		hdu = fits.PrimaryHDU(self.data)
	else:
		hdu = fits.PrimaryHDU(self.data.astype(np.float32))

	#Generate a header
	if hasattr(self,"cosmology") and (self.cosmology is not None):
		hdu.header["H0"] = (self.cosmology.H0.to(u.km/(u.s*u.Mpc)).value,"Hubble constant in km/s*Mpc")
		hdu.header["h"] = (self.cosmology.h,"Dimensionless Hubble constant")
		hdu.header["OMEGA_M"] = (self.cosmology.Om0,"Dark Matter density")
		hdu.header["OMEGA_L"] = (self.cosmology.Ode0,"Dark Energy density")
		hdu.header["W0"] = (self.cosmology.w0,"Dark Energy equation of state")
		hdu.header["WA"] = (self.cosmology.wa,"Dark Energy running equation of state")

	if hasattr(self,"redshift"):
		hdu.header["Z"] = (self.redshift,"Redshift of the background sources")
	
	hdu.header["ANGLE"] = (self.side_angle.to(u.deg).value,"Side angle in degrees")
	hdu.header["ITYPE"] = (self.__class__.__name__,"Image type") 

	#Save the image
	hdulist = fits.HDUList([hdu])
	hdulist.writeto(filename,clobber=True)