###############################################################
############CMB lensing reconstruction#########################
###############################################################
from __future__ import division,with_statement

import sys,os,glob

import numpy as np

from lenstools.image.convergence import ConvergenceMap,CMBTemperatureMap,PhiMap
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
		logdriver.info("CMB spectra type: {0}".format(settings.ps_type))
		logdriver.info("Reconstruction details: lmax={0} , sigmaN={1} , beam FWHM={2}".format(settings.lmax,settings.noise_level,settings.beam_fwhm))
		logdriver.info("Wiener filtering of reconstructed maps: {0}".format(settings.wiener))

	#Build noise data structure
	noise_keys = { "kind":"detector","sigmaN":settings.noise_level,"fwhm":settings.beam_fwhm }

	#Filtering
	if settings.wiener:
		filtering = "wiener"
	else:
		filtering = None

	############################
	#Match names of input files#
	############################

	#Find files, sort
	glob_pattern = os.path.join(maps_input.storage,settings.input_filename) 
	input_filenames = glob.glob(glob_pattern)
	input_filenames.sort()

	#Save information to info file for book-keeping
	if (pool is None) or (pool.is_master()):
		with open(os.path.join(maps_output.storage,"info.txt"),"w") as fp:
			fp.write("Inputs: {0}\n".format(glob_pattern))
			fp.write("Number of maps: {0}\n".format(len(input_filenames)))
			fp.write("Noise level: {0}\n".format(settings.noise_level))
			fp.write("Beam FWHM: {0}\n".format(settings.beam_fwhm))
			fp.write("lmax: {0}\n".format(settings.lmax))
			fp.write("Wiener filter: {0}\n".format(settings.wiener))

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
	realizations_per_task = len(input_filenames)//Ntask
	
	if pool is not None:
		first_realization = realizations_per_task*pool.rank
		last_realization = realizations_per_task*(pool.rank+1)
	else:
		first_realization = 0
		last_realization = len(input_filenames)

	#Log
	if pool is not None:
		logdriver.info("Task {0} processes inputs from {1} to {2}".format(pool.rank+1,first_realization+1,last_realization))

	#Conquer
	for n,r in enumerate(range(first_realization,last_realization)):

		#Set random seed
		np.random.seed(r)

		#Log
		input_filename = input_filenames[r]
		logdriver.debug("Loading input kappa image: {0}".format(input_filename))

		#Load input
		input_kappa = ConvergenceMap.load(input_filename)

		#Generate unlensed CMB sky
		if settings.estimator=="TT":
			logdriver.debug("Generating unlensed CMB sky...")
			cmb_unl = CMBTemperatureMap.from_power(angle=input_kappa.side_angle,npixel=input_kappa.data.shape[0],powerTT=unlensed_ps_filename,callback=settings.ps_type,lmax=settings.lmax)
		else:
			raise NotImplementedError("Quadratic estimator {0} not implemented!".format(settings.estimator))

		#Lens the CMB map, apply detector noise
		logdriver.debug("Lensing the CMB sky with input potential...")
		cmb_len = cmb_unl.lens(input_kappa)

		#Save intermediate step
		if settings.save_intermediate:
			logdriver.info("Saving intermediate steps: unlensed CMB sky, lensed CMB sky")
			cmb_unl.save(os.path.join(maps_output.storage,"cmb_unl_{0}".format(os.path.basename(input_filename))))
			cmb_len.save(os.path.join(maps_output.storage,"cmb_lensed_{0}".format(os.path.basename(input_filename))))

		logdriver.debug("Adding detector effects (white noise + beam deconvolution)...")
		cmb_len.addDetectorEffects(sigmaN=settings.noise_level,fwhm=settings.beam_fwhm)

		#Save intermediate step
		if settings.save_intermediate:
			logdriver.info("Saving intermediate steps: lensed CMB sky + detector effects")
			cmb_len.save(os.path.join(maps_output.storage,"cmb_lensed_detector_{0}".format(os.path.basename(input_filename))))

		#Apply quadratic estimator to reconstruct the lensing potential
		logdriver.debug("Estimating lensing {0} with {1} quadratic estimator...".format(settings.output_type,settings.estimator))

		if settings.output_type=="kappa":
			output_img = cmb_len.estimateKappaQuad(powerTT=lensed_ps_filename,callback=settings.ps_type,noise_keys=noise_keys,lmax=settings.lmax,filtering=filtering)

		elif settings.output_type=="phi":
			output_img = cmb_len.estimatePhiQuad(powerTT=lensed_ps_filename,callback=settings.ps_type,noise_keys=noise_keys,lmax=settings.lmax,filtering=filtering)

		elif settings.output_type=="phifft":
			phifft = cmb_len.estimatePhiFFTQuad(powerTT=lensed_ps_filename,callback=settings.ps_type,noise_keys=noise_keys,lmax=settings.lmax,filtering=filtering)
			output_img = PhiMap(phifft,angle=input_kappa.side_angle)
			
		else:
			raise NotImplementedError("output_type must be one between (kappa,phi,phifft)")

		#Save the result
		savename = os.path.join(maps_output.storage,settings.output_fname.format(os.path.basename(input_filename)))
		
		logdriver.info("Saving the reconstruction to: {0}".format(savename))
		output_img.save(savename)
		logdriver.debug("Saved the reconstruction to: {0}".format(savename))

		#Safety barrier
		if pool is not None:
			pool.comm.Barrier()

		#Log progress, memory usage
		peak_memory_task,peak_memory_all = peakMemory(),peakMemoryAll(pool)
		if (pool is None) or (pool.is_master()):
			logstderr.info("Progress: {0:.2f}%, peak memory usage: {1:.3f} (task), {2[0]:.3f} (all {2[1]} tasks)".format(100*(n+1)/realizations_per_task,peak_memory_task,peak_memory_all))

	#Complete
	if (pool is None) or (pool.is_master()):
		logdriver.info("Reconstructions complete, exiting.")
