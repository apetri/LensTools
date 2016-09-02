"""

.. module:: limber 
	:platform: Unix
	:synopsis: This module implements the tools to compute the convergence power spectrum from the 3D matter power spectrum using the Limber approximation


.. moduleauthor:: Andrea Petri <apetri@phys.columbia.edu>


"""

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

		self.kappa = self.kappa * kappa_units
		self.power = self.power * power_units

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