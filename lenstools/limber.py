import numpy as np
from scipy import interpolate,integrate

from astropy.constants import c

##################################################
#############Limber Integrator class##############
##################################################

class LimberIntegrator:
	
	"""
	A 3D power spectrum integrator that will compute the convergence power spectrum 
	using the Limber approximation.

	:param cosmoModel:
		One of astropy.cosmology objects

	:param lValues:
		Desired multipole values for the convergence power spectrum (numpy array)

	"""

	def __init__(self,cosmoModel,lValues):

		self.cosmoModel = cosmoModel
		self.lValues = lValues

	def computeConvergence(self,z,matterPower=None,powFileRoot=None):
		"""
		Computes the convergence power spectrum with the Limber integral of the 3D matter power spectrum.

		:param z:
			redshift bins at which the matter power spectrum is calculated

		:param matterPower:
			values of the matter power spectrum at corresponding z (first column must be k, the rest P(k,z), one for each z)

		:param powFileRoot:
			common root name of files in which the 3d power spectrum is stored; if None it is assumed that all the
			information is already loaded in the matterPower array. Throws and exception if both are None
		"""

		#Check validity of imput arguments
		if(matterPower==None and powFileRoot==None):
			raise ValueError("matterPower and powFileRoot cannot be both None!!")

		#Power spectrum normalization
		normalization = (9.0/4)*(self.cosmoModel.Om0)**2*(self.cosmoModel.H0.value/c.to("km/s").value)**4

		#l bins and convergence power spectrum
		l = self.lValues

		#Compute comoving distances and integral kernel
		chi = cosmoModel.comoving_distance(z)
		chi0 = chi[len(chi)-1]
		kernel = (1.0 - chi/chi0)**2

		#Load information about kappa (wavenumber) and P(kappa) (matter power spectrum)
		if(powFileRoot!=None):

			#Load matter power spectra from camb output files#
			#See how many k's are stored
			kappa,try_power = (np.loadtxt(powFileRoot + '0.dat')).transpose()
			num_kappa = len(kappa)

			#Load power spectrum
			power_spectrum = np.zeros([num_kappa,num_redshifts])

			for i in range(len(z)):
				try_power = np.loadtxt(powFileRoot + ('%d.dat'%(int(z[i]*100))))

				#Normalize power spectrum correctly
				power_spectrum[:,i] = try_power[:,1] / (cosmoModel.h**3)

		else:

			#Fill in values from the matterPower array
			kappa = matterPower[:,0]
			power_spectrum = matterPower[:,1:]

		#############################################################
		#Compute the integral for lensing convergence power spectrum#
		#############################################################
		
		power_interpolation = interpolate.interp1d(kappa,power_spectrum,axis=0)

		power_integrand = np.zeros((len(l),len(z)))
		lchi = (l/chi[np.newaxis,:]).reshape(len(l)*len(z))

		power_integrand = power_interpolation(lchi).reshape(len(l),len(z),len(z)).diagonal(axis1=1,axis2=2)
		full_integrand = kernel[np.newaxis,:] * (1.0 + z[np.newaxis,:])**2 * power_integrand
	
		#Finally compute the integral
		C = integrate.simps(full_integrand,chi,axis=1) * normalization

		#Return the final result
		return C


