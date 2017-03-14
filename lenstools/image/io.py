import numpy as np
from astropy.io import fits
import astropy.units as u

#####################################
#############FITS####################
#####################################

#Load scalar FITS image
def loadFITS(cls,filename):

	with fits.open(filename) as hdu:
		data = hdu[0].data

		if "ANGLE" in hdu[0].header:
			angle = hdu[0].header["ANGLE"] * u.deg
			return cls(data,angle)

		elif "SIDE" in hdu[0].header:
			side = hdu[0].header["SIDE"] * u.Mpc
			return cls(data,side)

		else:
			raise AttributeError("No ANGLE or SIDE keyword in header!")

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
	
	if self.side_angle.unit.physical_type=="angle":
		hdu.header["ANGLE"] = (self.side_angle.to(u.deg).value,"Side angle in degrees")
	else:
		hdu.header["SIDE"] = (self.side_angle.to(u.Mpc).value,"Side in Mpc")

	if hasattr(self,"unit"):
		hdu.header["UNIT"] = (self.unit.to_string(),"Pixel value unit")

	hdu.header["ITYPE"] = (self.__class__.__name__,"Image type") 

	#Save the image
	hdulist = fits.HDUList([hdu])
	hdulist.writeto(filename,overwrite=True)

#####################################
#############NUMPY###################
#####################################

#Write npz file
def saveNPZ(self,filename):
	np.savez(filename,self.data,angle_deg=self.side_angle.to(u.deg).value)

#Load npz file
def loadNPZ(cls,filename):
	data = np.load(filename)
	return cls(data["arr_0"],angle=data["angle_deg"]*u.deg)

