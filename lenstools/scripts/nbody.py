###################################################################################
################Measure statistics out of N-body simulation snapshots##############
###################################################################################

import sys,os
import logging

from lenstools.utils import MPIWhirlPool

from lenstools.simulations.nbody import NbodySnapshot

from lenstools.pipeline.simulation import SimulationBatch
from lenstools.pipeline.settings import PowerSpectrumSettings

import numpy as np

################################################
###########Loggers##############################
################################################

console = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter("%(asctime)s:%(name)-12s:%(levelname)-4s: %(message)s",datefmt='%m-%d %H:%M')
console.setFormatter(formatter)

logdriver = logging.getLogger("lenstools.driver")
logdriver.addHandler(console)
logdriver.propagate = False

################################################################
################Snapshot power spectrum#########################
################################################################

def powerSpectrum(pool,batch,settings,id,fmt):

	#Safety type check
	assert isinstance(pool,MPIWhirlPool) or (pool is None)
	assert isinstance(batch,SimulationBatch)
	assert isinstance(fmt(),NbodySnapshot)
	assert isinstance(settings,PowerSpectrumSettings)

	#Split the id into the model,collection and realization parts
	cosmo_id,geometry_id = id.split("|")

	#Get a handle on the simulation model
	model = batch.getModel(cosmo_id)

	#Scale the box size to the correct units
	nside,box_size = geometry_id.split("b")
	box_size = float(box_size)*model.Mpc_over_h

	#Get the handle on the collection
	collection = model.getCollection(box_size,nside)

	#Log the power spectrum settings to the user
	if (pool is None) or (pool.is_master()):

		logdriver.info("Measuring power spectrum for Ensemble {0}".format(settings.ensemble_name))
		logdriver.info("The Ensemble will be built with the following N-body realizations: {0}".format("-".join([ str(n) for n in settings.nbody_realizations ])))
		logdriver.info("First snapshot: {0}".format(settings.first_snapshot))
		logdriver.info("Last snapshot: {0}".format(settings.last_snapshot))
		logdriver.info("Minimum wavenumber: {0}".format(settings.kmin.to(model.Mpc_over_h**-1)))
		logdriver.info("Maximum wavenumber: {0}".format(settings.kmax.to(model.Mpc_over_h**-1)))
		logdriver.info("Bin size: {0}".format(((settings.kmax-settings.kmin)/settings.num_k_bins).to(model.Mpc_over_h**-1)))
		logdriver.info("FFT grid size: {0}".format(settings.fft_grid_size))
		logdriver.info("Number of bins: {0}".format(settings.num_k_bins))

	#Construct the array of bin edges
	k_egdes  = np.linspace(settings.kmin,settings.kmax,settings.num_k_bins+1).to(model.Mpc_over_h**-1)

	#Cycle over snapshots
	for n in range(settings.first_snapshot,settings.last_snapshot+1):

		#Allocate memory for the power spectrum ensemble
		power_ensemble = np.zeros((len(settings.nbody_realizations),settings.num_k_bins)) * (model.Mpc_over_h**3)

		#Log to user
		if (pool is None) or (pool.is_master()):
			logdriver.info("Processing snapshot {0} of model {1}".format(n,id))
			logdriver.info("Allocated memory for power spectrum Ensemble {0}".format(power_ensemble.shape))

		for r,ic in enumerate(settings.nbody_realizations):

			#Log to user
			if (pool is None) or (pool.is_master()):
				logdriver.info("Processing N-body realization {0}".format(ic))
			
			#Instantiate the appropriate SimulationIC object
			realization = collection.getRealization(ic)

			#Open the snapshot, measure the power spectrum and close it
			if pool is not None:
				if realization.gadget_settings.NumFilesPerSnapshot!=pool.size+1:
					logdriver.error("The number of snapshots written in parallel {0} does not coincide with the number of MPI processes {1}!".format(realization.gadget_settings.NumFilesPerSnapshot,pool.size+1))
					sys.exit(1)

			snap = fmt.open(realization.snapshotPath(n,sub=None),pool=pool)
			k,power_ensemble[r],hits = snap.powerSpectrum(k_egdes,resolution=settings.fft_grid_size,return_num_modes=True)
			snap.close()

			#Safety barrier sync
			if pool is not None:
				pool.comm.Barrier() 

		#Save the bin edges and mode counts
		if n==settings.first_snapshot and (pool is None or pool.is_master()):

			savename = os.path.join(collection.home_subdir,settings.ensemble_name+"_k.npy")
			logdriver.info("Saving wavevectors ({0}) to {1}".format(k.unit.to_string(),savename))
			np.save(savename,k.value)

			savename = os.path.join(collection.home_subdir,settings.ensemble_name+"_num_modes.npy")
			logdriver.info("Saving number of modes to {0}".format(savename))
			np.save(savename,hits)

		#Save the ensemble
		if (pool is None) or (pool.is_master()):
			
			savename = os.path.join(collection.home_subdir,settings.ensemble_name+"_snap{0:03d}.npy".format(n))
			logdriver.info("Saving power spectrum Ensemble ({0}) to {1}".format(power_ensemble.unit.to_string(),savename))
			np.save(savename,power_ensemble.value)


		#Safety barrier sync
		if pool is not None:
			pool.comm.Barrier()


	#Completed
	if pool is None or pool.is_master():
		logdriver.info("DONE!!")

