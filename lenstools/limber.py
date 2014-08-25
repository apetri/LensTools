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
from astropy.units import Mpc,def_unit

##################################################
#############Limber Integrator class##############
##################################################

class LimberIntegrator(object):
	
	"""
	A 3D power spectrum integrator that will compute the convergence power spectrum  using the Limber approximation.

	:param cosmoModel: One of astropy.cosmology objects (WMAP9 cosmology is set by default)

	:type cosmoModel: astropy.cosmology 

	"""

	def __init__(self,cosmoModel=cosmology.WMAP9):

		assert isinstance(cosmoModel,cosmology.FLRW),"cosmoModel should be a valid astropy cosmology instance!"
		self.cosmoModel = cosmoModel

		#Define also the Mpc/h units for convenience
		self.Mpc_over_h = def_unit("Mpc/h",Mpc/self.cosmoModel.h)

	def load3DPowerSpectrum(self,loader,*args,**kwargs):

		"""
		Loads in the matter power spectrum from (pre-computed) external files; args and kwargs are passed to the loader

		:param loader: must return, in order, k, z, P(k,z)
		:type loader: function

		"""

		self.z,self.kappa,self.power = loader(*args,**kwargs)
		assert self.power.shape == self.kappa.shape + self.z.shape
		assert self.z[0] > 0.0,"first redshift must be >0!!"

		self.setUnits()

	def setUnits(self,kappa_units=None,power_units=None):

		"""
		Set the physical units for wavenumber and matter power spectrum, default for length is Mpc

		"""

		if kappa_units is None:
			kappa_units = self.Mpc_over_h**-1

		if power_units is None:
			power_units = self.Mpc_over_h**3

		assert (power_units*(kappa_units**3)).physical_type == "dimensionless"

		assert hasattr(self,"kappa")
		assert hasattr(self,"power")

		self.kappa *= kappa_units
		self.power *= power_units

	def convergencePowerSpectrum(self,l):
		"""
		Computes the convergence power spectrum with the Limber integral of the 3D matter power spectrum; this still assumes a single source redshift at z0 = max(z) 
	
		:param l: multipole moments at which to compute the convergence power spectrum

		:type l: array

		:returns: array -- the convergence power spectrum at the l values specified

		"""

		z = self.z

		#Power spectrum normalization
		normalization = (9.0/4)*(self.cosmoModel.Om0)**2*(self.cosmoModel.H0/c)**4

		#Compute comoving distances and integral kernel (modify in case there is an arbitrary galaxy distribution)
		chi = self.cosmoModel.comoving_distance(z)
		chi0 = chi[-1]
		kernel = (1.0 - chi/chi0)**2

		#############################################################
		#Compute the integral for lensing convergence power spectrum#
		#############################################################
		
		power_interpolation = interpolate.interp1d(self.kappa.to(chi.unit**-1),self.power,axis=0,bounds_error=False,fill_value=0.0)
		lchi = np.outer(l,1.0/chi).reshape(len(l)*len(z))

		power_integrand = power_interpolation(lchi).reshape(len(l),len(z),len(z)).diagonal(axis1=1,axis2=2) * self.power.unit
		full_integrand = kernel[np.newaxis,:] * (1.0 + z[np.newaxis,:])**2 * power_integrand
	
		#Finally compute the integral
		C = integrate.simps(full_integrand,chi,axis=1) * normalization * full_integrand.unit * chi.unit
		assert C.unit.physical_type == u"dimensionless"

		#Return the final result
		return C.decompose().value

	def writeCAMBSettings(self,l,z,powFileRoot="matterpower",transfer_high_precision=False,transfer_k_per_logint=0,transfer_interp_matterpower=True):
		"""

		Outputs a StringIO object that will contain the redshift settings of the CAMB parameter file that will needed in order for CAMB to produce the linear or non linear matter power spectra that will then be integrated by the computeConvergence() method

		:param l: multipole moments at which to compute the convergence power spectrum
		:type l: array

		:param z: redshift bins at which the matter power spectrum is calculated (assumed to be a monotonic array with more than 1 element)

		:type z: array
		
		:param powFileRoot: root of the filename that you want to give to the CAMB power spectrum outputs

		:type powFileRoot: str.

		
		:param transfer_high_precision: read CAMB documentation (this sets the precision of the calculated transfer function)

		:type transfer_high_precision: bool.

		
		:param transfer_k_per_logint: read CAMB documentation (this sets the k wavenumber binning)

		:type transfer_k_per_logint: int.
		
		:param transfer_interp_matterpower: read CAMB documentation (this sets how the matter power is interpolated between different k's)

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






