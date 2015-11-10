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

###############################
######Multiquadric kernel######
###############################

def multiquadric(x,e):
	return np.sqrt(1.+x/e)

##################################
##Log of the gaussian likelihood##
##################################

def lnprobgauss(p,emulator,data,icov):
	diff = data - emulator.predict(p,raw=True)
	return -0.5*diff.dot(icov).dot(diff)

###############################
######emcee sampler############
###############################

def emcee_sampler(emulator,observed_feature,features_covariance,nwalkers=16,nburn=100,nchain=1000,pool=None):

	#Train the emulator with fast interpolation
	emulator.train(method=multiquadric)

	#Parameter space to sample
	parameters = emulator["parameters"].values
	pmin,pmax = parameters.min(0),parameters.max(0)

	#Feature name
	feature_name = emulator.feature_names[0]

	#Compute the inverse covariance
	icov = np.linalg.inv(features_covariance)

	#Initialize the walkers positions
	ndim = len(pmin)
	p0 = pmin + np.random.uniform(size=(nwalkers,ndim))*(pmax-pmin)

	#Initialize the sampler
	sampler = emcee.EnsembleSampler(nwalkers,ndim,lnprobgauss,args=[emulator,observed_feature,icov],pool=pool)

	#Burn-in
	logdriver.info("Running emcee burn-in: feature name={0}, feature dimension={1}, parameter dimension={2}, steps={3}".format(feature_name,len(observed_feature),ndim,nburn))
	pos,prob,state = sampler.run_mcmc(p0,nburn)
	sampler.reset()

	#Sampling
	logdriver.info("Running emcee MCMC chain: feature name={0}, feature dimension={1}, parameter dimension={2}, steps={3}".format(feature_name,len(observed_feature),ndim,nchain))
	sampler.run_mcmc(pos,nchain,rstate0=state)

	#Return sampled parameters
	return Ensemble(sampler.flatchain,columns=emulator.parameter_names)