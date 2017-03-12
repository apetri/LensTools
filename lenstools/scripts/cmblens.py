###############################################################
############CMB lensing reconstruction#########################
###############################################################
from __future__ import division,with_statement

import sys,os,glob

import numpy as np

from lenstools.simulations.logs import logdriver,logstderr,peakMemory,peakMemoryAll
from lenstools.pipeline.simulation import SimulationBatch
from lenstools.pipeline.settings import CMBReconstructionSettings

from lenstools.utils.mpi import MPIWhirlPool

##########################################
#Reconstruction with quadratic estimators#
##########################################

def reconstructQuad(pool,batch,settings,batch_id):

	#Safety check
	assert isinstance(pool,MPIWhirlPool) or (pool is None)
	assert isinstance(batch,SimulationBatch)
	assert isinstance(settings,CMBReconstructionSettings)

	#Parse cosmo_id,geometry
	cosmo_id,geometry_id = batch_id.split("|")

	#Get a hold of cosmological model, collection
	model = batch.getModel(cosmo_id)
	collection = model.getCollection(geometry_id)

	#Get a hold of input/output map sets
	maps_input = collection.getMapSet(settings.input_set)
	maps_output = collection.getMapSet(settings.output_set)

	#Log important settings to user
	if (pool is None) or (pool.is_master()):
		logdriver.info("Input map set with WL maps: {0}".format(settings.input_set))
		logdriver.info("Input file names match pattern: {0}".format(settings.input_filename))
		logdriver.info("Reconstructed lensing potential images will be saved to map set: {0}".format(settings.output_set))
		logdriver.info("Quadratic estimator type: {0} --> {1}".format(settings.estimator,settings.output_type))

	#Path to tabulated CMB power spectra
	unlensed_ps_filename = os.path.join(collection.home,settings.unlensed_ps_filename)
	lensed_ps_filename = os.path.join(collection.home,settings.lensed_ps_filename)

	#Log to user
	if (pool is None) or (pool.is_master()):
		logdriver.info("Unlensed CMB power spectra will be read from: {0}".format(unlensed_ps_filename))
		logdriver.info("Lensed CMB power spectra will be read from: {0}".format(lensed_ps_filename))
		logdriver.info("Reconstruction details: lmax={0} , sigmaN={1} , beam FWHM={2}".format(settings.lmax,settings.noise_level,settings.beam_fwhm))
		logdriver.info("Wiener filtering of reconstructed maps: {0}".format(settings.wiener))

	############################
	#Match names of input files#
	############################

	#Find files, sort
	input_filenames = glob.glob(os.path.join(maps_input.storage,settings.input_filename))
	input_filenames.sort()

	#Divide uniformly between tasks
	if pool is None:
		Ntask = 1
	else:
		Ntask = pool.size+1

	#Log 
	if (pool is None) or (pool.is_master()):
		logdriver.info("Found {0} inputs matching pattern {1}, will be divided between {2} tasks".format(len(input_filenames),settings.input_filename,Ntask))

	#Perfect load balance enforced
	if len(input_filenames)%Ntask:
		if (pool is None) or (pool.is_master()):
			logdriver.error("Perfect load balance enforced: Ntask={0} must divide the number of inputs ({1})".format(Ntask,len(input_filenames)))
		sys.exit(1)

	#Divide
	if pool is not None:
		realizations_per_task = len(input_filenames)//Ntask
		first_realization = realizations_per_task*pool.rank
		last_realization = realizations_per_task*(pool.rank+1)
	else:
		first_realization = 0
		last_realization = len(input_filenames)

	#Log
	if pool is not None:
		logdriver.info("Task {0} processes inputs from {1} to {2}".format(pool.rank+1,first_realization+1,last_realization))

	#Conquer
	for r in range(first_realization,last_realization):

		#Set random seed
		np.random.seed(r)

		#Log
		input_filename = input_filenames[r]
		logdriver.info("Processing {0}".format(input_filename))
