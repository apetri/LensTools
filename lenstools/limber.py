"""

.. module:: limber 
	:platform: Unix
	:synopsis: This module implements the tools to compute the convergence power spectrum from the 3D matter power spectrum using the Limber approximation


.. moduleauthor:: Andrea Petri <apetri@phys.columbia.edu>


"""

import StringIO

import numpy as np
from scipy import interpolate,integrate

from astropy import cosmology
from astropy.constants import c

##################################################
#############Limber Integrator class##############
##################################################

class LimberIntegrator(object):
	
	"""
	A 3D power spectrum integrator that will compute the convergence power spectrum  using the Limber approximation. The units for quantities with dimensions of length are assumed to be Mpc

	
	:param lValues: Desired multipole values for the convergence power spectrum

	:type lValues: array

	:param cosmoModel: One of astropy.cosmology objects (WMAP9 cosmology is set by default)

	:type cosmoModel: astropy.cosmology 

	"""

	def __init__(self,lValues,cosmoModel=cosmology.WMAP9):

		self.lValues = lValues
		self.cosmoModel = cosmoModel

	def computeConvergence(self,z,matterPower=None,powFileRoot=None,extension=".dat"):
		"""
		Computes the convergence power spectrum with the Limber integral of the 3D matter power spectrum;
		this still assumes a single source redshift at z0 = max(z) 

		
		:param z:
			redshift bins at which the matter power spectrum is calculated; the z array must be sorted in ascending order
			and z[0] must be >0 even if small, otherwise this throws an exception

		:type z: array

		
		:param matterPower:
			values of the matter power spectrum at corresponding z (first column must be k, the rest P(k,z), one for each z)

		:type matterPower: ndarray

		
		:param powFileRoot:
			common root name of files in which the 3d power spectrum is stored; if None it is assumed that all the
			information is already loaded in the matterPower array. Throws and exception if both are None

		:type powFileRoot: str.

		
		:param extension:
			extension of text files with 3d power spectrum, default is .dat

		:type extension: str.

		:returns: array -- the convergence power spectrum at the l Values specified in the constructor

		:raises: ValueError
		

		"""

		#Check validity of redshift values
		if(z[0]<=0.0):
			raise ValueError("first redshift must be >0!!")

		#Check validity of imput arguments
		if(matterPower==None and powFileRoot==None):
			raise ValueError("matterPower and powFileRoot cannot be both None!!")

		#Power spectrum normalization
		normalization = (9.0/4)*(self.cosmoModel.Om0)**2*(self.cosmoModel.H0.value/c.to("km/s").value)**4

		#l bins and convergence power spectrum
		l = self.lValues

		#Compute comoving distances and integral kernel
		chi = self.cosmoModel.comoving_distance(z)
		chi0 = chi[len(chi)-1]
		kernel = (1.0 - chi/chi0)**2

		#Load information about kappa (wavenumber) and P(kappa) (matter power spectrum)
		if(powFileRoot!=None):

			#Load matter power spectra from camb output files#
			#See how many k's are stored
			kappa,try_power = (np.loadtxt(powFileRoot + int(z[0]*100) + extension)).transpose()

			#Load power spectrum
			power_spectrum = np.zeros([len(kappa),len(z)])

			for i in range(len(z)):
				try_power = np.loadtxt(powFileRoot + ('%d%s'%(int(z[i]*100,extension))))

				#Normalize power spectrum correctly
				power_spectrum[:,i] = try_power[:,1] / (self.cosmoModel.h**3)

		else:

			#Fill in values from the matterPower array
			kappa = matterPower[:,0]
			power_spectrum = matterPower[:,1:]

		#############################################################
		#Compute the integral for lensing convergence power spectrum#
		#############################################################
		
		power_interpolation = interpolate.interp1d(kappa*self.cosmoModel.h,power_spectrum,axis=0,bounds_error=False,fill_value=0.0)

		power_integrand = np.zeros((len(l),len(z)))
		lchi = np.outer(l,1.0/chi).reshape(len(l)*len(z))

		power_integrand = power_interpolation(lchi).reshape(len(l),len(z),len(z)).diagonal(axis1=1,axis2=2)
		full_integrand = kernel[np.newaxis,:] * (1.0 + z[np.newaxis,:])**2 * power_integrand
	
		#Finally compute the integral
		C = integrate.simps(full_integrand,chi,axis=1) * normalization

		#Return the final result
		return C

	def writeCAMBSettings(self,z,powFileRoot="matterpower",transfer_high_precision=False,transfer_k_per_logint=0,transfer_interp_matterpower=True):
		"""
		
		Outputs a StringIO object that will contain the redshift settings of the CAMB parameter file that will needed
		in order for CAMB to produce the linear or non linear matter power spectra that will then be integrated by 
		the computeConvergence() method

		
		:param z:
			redshift bins at which the matter power spectrum is calculated (assumed to be a monotonic array with more than 1 element)

		:type z: array

		
		:param powFileRoot:
			root of the filename that you want to give to the CAMB power spectrum outputs

		:type powFileRoot: str.

		
		:param transfer_high_precision:
			read CAMB documentation (this sets the precision of the calculated transfer function)

		:type transfer_high_precision: bool.

		
		:param transfer_k_per_logint:
			read CAMB documentation (this sets the k wavenumber binning)

		:type transfer_k_per_logint: int.
		
		:param transfer_interp_matterpower:
			read CAMB documentation (this sets how the matter power is interpolated between different k's)

		:type transfer_interp_matterpower: bool.

		:returns: str -- the portion of the CAMB parameter file relevant to the 3D matter power spectrum

		"""

		S = StringIO.StringIO()

		if(transfer_high_precision):
			S.write("""transfer_high_precision = T\n""")
		else:
			S.write("""transfer_high_precision = F\n""")

		#Compute maximum k that needs to be calculated (tune manually in the parameter file if this is too high)
		if(z[0]==0.0):
			kmax = self.lValues.max()/self.cosmoModel.comoving_distance(z[1]).value
		else:
			kmax = self.lValues.max()/self.cosmoModel.comoving_distance(z[0]).value

		S.write("""transfer_kmax = %.3f
transfer_k_per_logint = %d
transfer_num_redshifts  = %d

"""%(kmax,transfer_k_per_logint,len(z)))

		#Other settings
		if(transfer_interp_matterpower):
			S.write("""transfer_interp_matterpower = T\n\n""")
		else:
			S.write("""transfer_interp_matterpower = F\n\n""")

		#Write desired output redshifts
		for i in range(len(z)):
			S.write("""transfer_redshift(%d) = %.3f\n"""%(i+1,z[i]))
			S.write("""transfer_filename(%d) = transfer_out_%d.dat\n"""%(i+1,int(z[i]*100)))

		#Write desired output filenames
		S.write("""\n#Matter power spectrum output against k/h in units of h^{-3} Mpc^3\n\n""")

		for i in range(len(z)):
			S.write("""transfer_matterpower(%d) = %s%d.dat\n"""%(i+1,powFileRoot,int(z[i]*100)))

		#Return the StringIO object
		S.seek(0)
		return S.read()






