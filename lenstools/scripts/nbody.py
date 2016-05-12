###################################################################################
################Measure statistics out of N-body simulation snapshots##############
###################################################################################

import sys,os
import logging

from distutils import config
from ConfigParser import NoOptionError

from lenstools.utils import MPIWhirlPool

from lenstools.simulations.nbody import NbodySnapshot
from lenstools.simulations.gadget2 import Gadget2SnapshotDE
from lenstools.simulations.logs import logdriver

from lenstools.pipeline.simulation import SimulationBatch

import numpy as np
import astropy.units as u

#Orchestra director of the execution
def powerSpectrumExecution():

	script_to_execute = matterPowerSpectrum
	settings_handler = PowerSpectrumSettings
	kwargs = {"fmt":Gadget2SnapshotDE}

	return script_to_execute,settings_handler,kwargs

################################################################
################Snapshot power spectrum#########################
################################################################

def matterPowerSpectrum(pool,batch,settings,node_id,**kwargs):

	assert "fmt" in kwargs.keys()
	fmt = kwargs["fmt"]

	#Safety type check
	assert isinstance(pool,MPIWhirlPool) or (pool is None)
	assert isinstance(batch,SimulationBatch)
	assert isinstance(fmt(),NbodySnapshot)
	assert isinstance(settings,PowerSpectrumSettings)

	#Split the id into the model,collection and realization parts
	cosmo_id,geometry_id = node_id.split("|")

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

		#Create dedicated ensemble directory
		ensemble_dir = os.path.join(collection.home_subdir,settings.ensemble_name)
		if not os.path.isdir(ensemble_dir):
			os.mkdir(ensemble_dir) 

	#Construct the array of bin edges
	k_egdes  = np.linspace(settings.kmin,settings.kmax,settings.num_k_bins+1).to(model.Mpc_over_h**-1)

	#Placeholder for the density MPI communications
	density_placeholder = np.empty((settings.fft_grid_size,)*3,dtype=np.float32)
	if pool is not None:
		pool.openWindow(density_placeholder)

		if pool.is_master():
			logdriver.debug("Opened density window of type {0}".format(pool._window_type))

	#Cycle over snapshots
	for n in range(settings.first_snapshot,settings.last_snapshot+1):

		#Allocate memory for the power spectrum ensemble
		power_ensemble = np.zeros((len(settings.nbody_realizations),settings.num_k_bins)) * (model.Mpc_over_h**3)

		#Log to user
		if (pool is None) or (pool.is_master()):
			logdriver.info("Processing snapshot {0} of model {1}".format(n,node_id))
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
			k,power_ensemble[r],hits = snap.powerSpectrum(k_egdes,resolution=settings.fft_grid_size,return_num_modes=True,density_placeholder=density_placeholder)
			snap.close()

			#Safety barrier sync
			if pool is not None:
				pool.comm.Barrier() 

		#Save the bin edges and mode counts
		if n==settings.first_snapshot and (pool is None or pool.is_master()):

			savename = os.path.join(collection.home_subdir,settings.ensemble_name,settings.ensemble_name+"_k.npy")
			logdriver.info("Saving wavevectors ({0}) to {1}".format(k.unit.to_string(),savename))
			np.save(savename,k.value)

			savename = os.path.join(collection.home_subdir,settings.ensemble_name,settings.ensemble_name+"_num_modes.npy")
			logdriver.info("Saving number of modes to {0}".format(savename))
			np.save(savename,hits)

		#Save the ensemble
		if (pool is None) or (pool.is_master()):
			
			savename = os.path.join(collection.home_subdir,settings.ensemble_name,settings.ensemble_name+"_snap{0:03d}.npy".format(n))
			logdriver.info("Saving power spectrum Ensemble ({0}) to {1}".format(power_ensemble.unit.to_string(),savename))
			np.save(savename,power_ensemble.value)


		#Safety barrier sync
		if pool is not None:
			pool.comm.Barrier()

	###########
	#Completed#
	###########

	#Close the RMA window
	if pool is not None:
		pool.comm.Barrier()
		pool.closeWindow()
		
		if pool.is_master():
			logdriver.debug("Closed density window of type {0}".format(pool._window_type))


	if pool is None or pool.is_master():
		logdriver.info("DONE!!")



########################################################
###########PowerSpectrumSettings class##################
########################################################

class PowerSpectrumSettings(object):

	"""
	Class handler of N-Body simulation power spectrum measurement settings

	"""

	def __init__(self,**kwargs):

		#Tunable settings (resolution, etc...)
		self.ensemble_name = "gadget2_ps"
		self.nbody_realizations = [1]
		self.first_snapshot = 46
		self.last_snapshot = 58
		self.fft_grid_size = 256
		self.kmin = 0.003 * u.Mpc**-1
		self.kmax = 1.536 * u.Mpc**-1
		self.length_unit = u.Mpc
		self.num_k_bins = 50

		#Allow for kwargs override
		for key in kwargs.keys():
			setattr(self,key,kwargs[key])

	@classmethod
	def read(cls,config_file):

		#Read the options from the ini file
		options = config.ConfigParser()
		options.read([config_file])

		#Check that the config file has the appropriate section
		section = "PowerSpectrumSettings"
		assert options.has_section(section),"No {0} section in configuration file {1}".format(section,config_file)

		#Fill in the appropriate fields
		settings = cls()

		settings.ensemble_name = options.get(section,"ensemble_name")

		#Read in the nbody realizations that make up the ensemble
		settings.nbody_realizations = list()
		for r in options.get(section,"nbody_realizations").split(","): 
			
			try:
				l,h = r.split("-")
				settings.nbody_realizations.extend(range(int(l),int(h)+1))
			except ValueError:
				settings.nbody_realizations.extend([int(r)])
		
		settings.first_snapshot = options.getint(section,"first_snapshot")
		settings.last_snapshot = options.getint(section,"last_snapshot")
		
		settings.fft_grid_size = options.getint(section,"fft_grid_size")

		settings.length_unit = getattr(u,options.get(section,"length_unit"))
		settings.kmin = options.getfloat(section,"kmin") * settings.length_unit**-1
		settings.kmax = options.getfloat(section,"kmax") * settings.length_unit**-1
		
		settings.num_k_bins = options.getint(section,"num_k_bins")

		#Return to user
		return settings

