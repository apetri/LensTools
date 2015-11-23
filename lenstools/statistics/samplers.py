"""

.. module:: samplers
	:platform: Unix
	:synopsis: This module contains wrappers for some pre implemented parameter samplers


.. moduleauthor:: Andrea Petri <apetri@phys.columbia.edu>


"""

from __future__ import division,print_function,with_statement

import numpy as np
import emcee

from .ensemble import Ensemble
from ..simulations.logs import logdriver
from ..utils.algorithms import precision_bias_correction

###############################
######Multiquadric kernel######
###############################

def multiquadric(x,e):
	return np.sqrt(1.+x/e)

##################################
##Log of the gaussian likelihood##
##################################

def lnprobgauss(p,emulator,data,icov,pslice_values,sample_indices):

	if pslice_values is not None:
		
		pc = np.empty(len(pslice_values)+len(sample_indices))
		for i in pslice_values:
			pc[i] = pslice_values[i]
		for n,i in enumerate(sample_indices):
			pc[i] = p[n]

	else:
		pc =p

	diff = data - emulator.predict(pc,raw=True)
	return -0.5*diff.dot(icov).dot(diff)

###############################
######emcee sampler############
###############################

def emcee_sampler(emulator,observed_feature,features_covariance,correct=None,pslice=None,nwalkers=16,nburn=100,nchain=1000,pool=None):

	"""
	Parameter posterior sampling based on the MCMC algorithm implemented by the emcee package

	:param emulator: feature emulator
	:type emulator: :py:class:`Emulator`

	:param observed_feature: observed feature to condition the parameter estimation
	:type observed_feature: array

	:param features_covariance: covariance matrix of the features
	:type features_covariance: array

	:param correct: if not None, correct for the bias in the inverse covariance estimator assuming the covariance was estimated by 'correct' simulations
	:type correct: int.

	:param pslice: specify slices of the parameter space in which some parameters are keps as constants
	:type pslice: dict.

	:param nwalkers: number of chains
	:type nwalkers: int.

	:param nburn: length of the burn-in chain
	:type nburn: int.

	:param nchain: length of the MCMC chain
	:type nchain: int.

	:param pool: MPI Pool for parallelization of computations
	:type pool: MPIPool

	:returns: ensemble of samples from the posterior probability distribution
	:rtype: :py:class:`Ensemble`

	"""

	#Train the emulator with fast interpolation
	emulator.train(method=multiquadric)
	parameters = emulator["parameters"]
	parameters_to_sample = emulator.parameter_names

	#Parameter space to sample
	if pslice is None: 
		parameters = parameters.values
		pslice_values = None
		sample_indices = None
	else:
		pslice_values = dict((emulator.parameter_names.index(p),pslice[p]) for p in pslice)
		sample_indices = [ emulator.parameter_names.index(p) for p in emulator.parameter_names if p not in pslice ]
		assert len(pslice_values)+len(sample_indices)==len(emulator.parameter_names)
		parameters = parameters.values[:,sample_indices]

	#Extremes of the sampling space	
	pmin,pmax = parameters.min(0),parameters.max(0)

	#Feature name
	feature_name = emulator.feature_names[0]

	#Compute the inverse covariance
	icov = np.linalg.inv(features_covariance)
	if correct is not None:
		icov *= precision_bias_correction(correct,len(icov))

	#Initialize the walkers positions
	ndim = len(pmin)
	p0 = pmin + np.random.uniform(size=(nwalkers,ndim))*(pmax-pmin)

	#Initialize the sampler
	sampler = emcee.EnsembleSampler(nwalkers,ndim,lnprobgauss,args=[emulator,observed_feature,icov,pslice_values,sample_indices],pool=pool)

	#Burn-in
	logdriver.info("Running emcee burn-in: feature name={0}, feature dimension={1}, parameter dimension={2}, steps={3}".format(feature_name,len(observed_feature),ndim,nburn))
	pos,prob,state = sampler.run_mcmc(p0,nburn)
	sampler.reset()

	#Sampling
	logdriver.info("Running emcee MCMC chain: feature name={0}, feature dimension={1}, parameter dimension={2}, steps={3}".format(feature_name,len(observed_feature),ndim,nchain))
	sampler.run_mcmc(pos,nchain,rstate0=state)

	#Return sampled parameters
	if pslice is None:
		return Ensemble(sampler.flatchain,columns=emulator.parameter_names)
	else:
		return Ensemble(sampler.flatchain,columns=filter(lambda p:p not in pslice,emulator.parameter_names))