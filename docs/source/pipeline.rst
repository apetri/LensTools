The LensTools Weak Lensing Simulation pipeline
**********************************************

This document is a hands--on tutorial on how to deploy the lenstools weak lensing simulation pipeline on a computer cluster. This will enable you to run your own weak lensing simulations and produce simulated weak lensing fields (shear and convergence) starting from a set of cosmological parameters. 

Pipeline workflow
-----------------

:ref:`flow` can be seen in the figure below. The starting point is to select a set of cosmological parameters that are supposed to describe your mock universe (i.e. a combination of :math:`(H_0,\Omega_m,\Omega_b,\Omega_\Lambda,w_0,w_a,\sigma_8,n_s)`). The first step in the pipeline is running several independent :math:`N`--body simulations that will give an accurate picture of the dark matter density field of the universe at different redshifts. Once the :math:`N`--body boxes are produced, the multi--lens--plane algorithm is used to trace the deflection of light rays as they go through the simulated density boxes. This will allow to compute background galaxies shape distortions and/or weak lensing quantities such as the shear :math:`\gamma_{1,2}` and the convergence :math:`\kappa` (and, with a limited accuracy, the rotation angle :math:`\omega`). 

.. _flow:

.. figure:: flow.png

	A scheme of the pipeline workflow 

In the remainder of the document a detailed scheme of the pipeline is illustrated.

The directory tree
------------------

lenstools provides a set of routines for managing the simulations directory tree, which is crucial for organizing the produced files in a sensible way. The first step in the process is creating a new batch of simulations. The way you accomplish this is through a SimulationBatch object. 

::
	
	>>> from lenstools.pipeline import SimulationBatch
	>>> from lenstools.pipeline.settings import EnvironmentSettings

You will need to choose where you want to store your files: in each simulation batch there are two distinct locations files will be saved in. The "home" location is reserved for small files such as code parameter files, tabulated power spectra and other book--keeping necessary files. The "storage" location is used to store large production files, such as :math:`N`--body simulation boxes, lensing planes and weak lensing maps. These locations need to be specified upon the batch creation

::

	environment = EnvironmentSettings(home="SimTest/Home",storage="SimTest/Storage")
	batch = SimulationBatch(environment)

You will need to specify the home and storage paths only once throughout the execution of the pipeline, lenstools will do the rest! If you want to build a git repository on top of your simulation batch, you will have to install `GitPython <https://gitpython.readthedocs.org>`_ and initiate the simulation batch as follows

::
	
	>>> from lenstools.pipeline.remote import LocalGit
	>>> batch = SimulationBatch(environment,syshandler=LocalGit())

Cosmological parameters
~~~~~~~~~~~~~~~~~~~~~~~

We first need to specify the cosmological model that will regulate the physics of the universe expansion and evolution of the density perturbations. We do this through the cosmology module of astropy (slightly modified to allow the specifications of parameters like :math:`n_s` and :math:`\sigma_8`). 

::

	>>> from lenstools.pipeline.simulation import CosmoDefault
	>>> cosmology = CosmoDefault(Om0=0.3,Ode0=0.7)

The cosmology object will be initialized with :math:`(\Omega_m,\Omega_\Lambda)=(0.3,0.7)` and all the other parameters set to their default values

::

	>>> cosmology
	Nicaea(H0=72 km / (Mpc s), Om0=0.3, Ode0=0.7, sigma8=0.8, ns=0.96, w0=-1, wa=0, Tcmb0=2.725 K, Neff=3.04, m_nu=[ 0.  0.  0.] eV, Ob0=0.046) 

Now we create a new simulation model that corresponds to the "cosmology" just specified, through our "batch" handler created before

:: 

	>>> model = batch.newModel(cosmology,parameters=["Om","Ol"])
	
	[+] SimTest/Home/Om0.300_Ol0.700 created
	[+] SimTest/Storage/Om0.300_Ol0.700 created

The argument "parameters" specifies which cosmological parameters you want to keep track of in your model; this is useful, for example, when you want to simulate different combinations of these parameters while keeping the other fixed to their default values. Note that lenstools informs you of the directories that are created on disk. You have access at any time to the models that are present in your simulation batch 

::

	>> batch.available
	[<Om=0.300 , Ol=0.700>]


Simulation resolution
~~~~~~~~~~~~~~~~~~~~~

It is now time to specify the resolution of the :math:`N`--body simulations that will be run to map the 3D density field of the universe. There are two numbers you need to set here, namely size of the box (that will fix the largest mode your simulations will be able to probe) and the number of particles on a side (that will fix the shortest mode). This command will create a collection of simulations with :math:`512^3` particles in a box of size 240.0 Mpc/h

::

	>>> collection = model.newCollection(box_size=240.0*model.Mpc_over_h,nside=512)
	
	[+] SimTest/Home/Om0.300_Ol0.700/512b240 created
	[+] SimTest/Storage/Om0.300_Ol0.700/512b240 created

Again, you will have access at any time to the collections that are present in your model 

::

	>>> model = batch.getModel("Om0.300_Ol0.700")
	>>> model.collections 
	
	[<Om=0.300 , Ol=0.700> | box=240.0 Mpc/h,nside=512]

Initial conditions
~~~~~~~~~~~~~~~~~~

Each simulation collection can have multiple realizations of the density field; these realizations share all the same statistical properties (i.e. the matter power spectrum), but have different spatial arrangements of the particles. This allows you to measure ensemble statistics such as means and covariances of various observables. Let's add three independent realizations of the density field to the "512b240" collection, with random seeds 1,22,333 (the random seed will be used by the initial condition generator to produce different density fields that share the same 3D power spectum)

::

	>>> for s in [1,22,333]:
		collection.newRealization(seed=s)

	[+] SimTest/Home/Om0.300_Ol0.700/ic1 created
	[+] SimTest/Storage/Om0.300_Ol0.700/ic1 created
	[+] SimTest/Home/Om0.300_Ol0.700/ic2 created
	[+] SimTest/Storage/Om0.300_Ol0.700/ic2 created
	[+] SimTest/Home/Om0.300_Ol0.700/ic3 created
	[+] SimTest/Storage/Om0.300_Ol0.700/ic3 created 

At this point it should not be surprising that you can do this

::

	>>> collection.realizations

	[<Om=0.300 , Ol=0.700> | box=240.0 Mpc/h,nside=512 | ic=1,seed=1 | IC files on disk: 0 | Snapshot files on disk: 0,
 	<Om=0.300 , Ol=0.700> | box=240.0 Mpc/h,nside=512 | ic=2,seed=22 | IC files on disk: 0 | Snapshot files on disk: 0,
 	<Om=0.300 , Ol=0.700> | box=240.0 Mpc/h,nside=512 | ic=3,seed=333 | IC files on disk: 0 | Snapshot files on disk: 0]

Note that, at this step, we are only laying down the directory tree of the simulation batch, and you can see that there are neither IC files nor snapshot files saved on disk yet (this will be produced when we actually run the simulations, but this will be explained later in the tutorial). 


Lens planes
~~~~~~~~~~~

For each of the realizations in the collection, we have to create a set of lens planes, that will be necessary for the execution of the ray--tracing step via the multi--lens--plane algorithm. The settings for these lens plane set can be specified through a INI configuration file. Let's call this file "planes.ini"; it should have the following structure

::

	[PlaneSettings]

	directory_name = Planes
	override_with_local = False
	format = fits
	plane_resolution = 128
	first_snapshot = 0
	last_snapshot = 58
	cut_points = 10.71
	thickness = 3.57 
	length_unit = Mpc
	normals = 0,1,2

Once you specified the plane configuration file, you can go ahead and create a lens plane set for each of the :math:`N`--body realizations you created at the previous step

::

	>>> from lenstools.pipeline.settings import PlaneSettings
	>>> plane_settings = PlaneSettings.read("planes.ini")
	>>> for r in collection.realizations:
		r.newPlaneSet(plane_settings)

	[+] SimTest/Home/Om0.300_Ol0.700/ic1/Planes created
	[+] SimTest/Storage/Om0.300_Ol0.700/ic1/Planes created
	[+] SimTest/Home/Om0.300_Ol0.700/ic2/Planes created
	[+] SimTest/Storage/Om0.300_Ol0.700/ic2/Planes created
	[+] SimTest/Home/Om0.300_Ol0.700/ic3/Planes created
	[+] SimTest/Storage/Om0.300_Ol0.700/ic3/Planes created

To summarize what you just did, as usual you can type 

::

	>>> for r in collection.realizations:
		r.planesets

	[<Om=0.300 , Ol=0.700>  |  box=240.0 Mpc/h,nside=512  |  ic=1,seed=1  | Plane set: Planes , Plane files on disk: 0]
	[<Om=0.300 , Ol=0.700>  |  box=240.0 Mpc/h,nside=512  |  ic=2,seed=22  | Plane set: Planes , Plane files on disk: 0]
	[<Om=0.300 , Ol=0.700>  |  box=240.0 Mpc/h,nside=512  |  ic=3,seed=3333  | Plane set: Planes , Plane files on disk: 0]


Weak lensing fields
~~~~~~~~~~~~~~~~~~~

The last step in the pipeline is to run the multi--lens--plane algorithm through the sets of lens planes just created. This will compute all the ray deflections at each lens crossing and derive the corresponding weak lensing quantities. The ray tracing settings need to be specified in a INI configuration file, that for example we can call "lens.ini". The following configuration will allow you to create square weak lensing simulated maps assuming all the background sources have the same redshift 

::

	[MapSettings]

	directory_name = Maps
	override_with_local = False
	format = fits
	map_resolution = 128
	map_angle = 3.5
	angle_unit = deg
	source_redshift = 2.0

	#Random seed used to generate multiple map realizations
	seed = 0

	#Set of lens planes to be used during ray tracing
	plane_set = Planes

	#N-body simulation realizations that need to be mixed
	mix_nbody_realizations = 1,2,3
	mix_cut_points = 0,1,2
	mix_normals = 0,1,2
	lens_map_realizations = 4

	#Which lensing quantities do we need?
	convergence = True
	shear = True
	omega = True

Different random realizations of the same weak lensing field can be obtained drawing different combinations of the lens planes from different :math:`N`--body realizations (*mix_nbody_realizations*), different regions of the :math:`N`--body boxes (*mix_cut_points*) and different rotation of the boxes (*mix_normals*). We create the directories for the weak lensing map set as usual

::

	>>> from lenstools.pipeline.settings import MapSettings
	>>> map_settings = MapSettings.read("lens.ini")
	>>> map_set = collection.newMapSet(map_settings)

	[+] SimTest/Home/Om0.300_Ol0.700/Maps created
	[+] SimTest/Storage/Om0.300_Ol0.700/Maps created

And, of course, you can check what you just did 

::

	>>> collection.mapsets

	[<Om=0.300 , Ol=0.700> | box=240.0 Mpc/h,nside=512 | Map set: Maps | Map files on disk: 0 ]

Now that we layed down our directory tree in a logical and organized fashion, we can proceed with the deployment of the simulation codes. The outputs of these codes will be saved in the "storage" portion of the simulation batch. 

Pipeline deployment
-------------------

Matter power spectra (CAMB)
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Initial conditions (NGenIC)
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Gravitational evolution (Gadget2)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Lens planes
~~~~~~~~~~~

Weak lensing fields :math:`\gamma,\kappa,\omega`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Post processing
---------------

Example: measure the 3D matter power spectrum
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Create a file "matter_power_spectrum.py" 

::

	###################################################################################
	################Measure statistics out of N-body simulation snapshots##############
	###################################################################################

	import sys,os
	import logging

	from distutils import config
	from ConfigParser import NoOptionError

	from lenstools.utils import MPIWhirlPool

	from lenstools.simulations.nbody import NbodySnapshot
	from lenstools.simulations.gadget2 import Gadget2Snapshot

	from lenstools.pipeline.simulation import SimulationBatch

	import numpy as np
	import astropy.units as u

	################################################
	###########Loggers##############################
	################################################

	console = logging.StreamHandler(sys.stdout)
	formatter = logging.Formatter("%(asctime)s:%(name)-12s:%(levelname)-4s: %(message)s",datefmt='%m-%d %H:%M')
	console.setFormatter(formatter)

	logdriver = logging.getLogger("lenstools.driver")
	logdriver.addHandler(console)
	logdriver.propagate = False

	#Orchestra director of the execution
	def powerSpectrumExecution():

		script_to_execute = matterPowerSpectrum
		settings_handler = PowerSpectrumSettings
		kwargs = {"fmt":Gadget2Snapshot}

		return script_to_execute,settings_handler,kwargs

	################################################################
	################Snapshot power spectrum#########################
	################################################################

	def matterPowerSpectrum(pool,batch,settings,id,**kwargs):

		assert "fmt" in kwargs.keys()
		fmt = kwargs["fmt"]

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

			#Create dedicated ensemble directory
			ensemble_dir = os.path.join(collection.home_subdir,settings.ensemble_name)
			if not os.path.isdir(ensemble_dir):
				os.mkdir(ensemble_dir) 

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


		#Completed
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


Create a INI configuration file "code_options.ini": 

::
	
	[PowerSpectrumSettings]

	ensemble_name = gadget2_ps
	nbody_realizations = 1,2-3
	first_snapshot = 0
	last_snapshot = 58
	fft_grid_size = 256
	kmin = 0.003 
	kmax = 1.536 
	length_unit = Mpc
	num_k_bins = 50

You deploy like this 

::

	lenstools.execute-mpi -e SimTest/environment.ini -c code_options.ini -m matter_power_spectrum -p powerSpectrumExecution "Om0.300_Ol0.700|512b240"












