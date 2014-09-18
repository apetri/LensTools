Weak lensing simulations
************************

At the current stage of Weak gravitational lensing research, large numerical simulations are required for analyzing observations; the LensTools python package provides an API to interact with some already existing simulated datasets (mainly convergence and shear maps for different cosmological models), such as 

1. The IGS1 simulations: this simulated dataset contains 1000 realizations of single redshift convergence and shear maps for six different cosmological parammeter combinations (a fiducial model and some variations of the quadruplet :math:`(\Omega_m,w,\sigma_8,n_s)`). The fiducial model is based on 45 independent N-body simulations and the variations are based on 5 independent N-body simulations (where :math:`N=512^3`)

2. The CFHTemu1 simulations: this simulated dataset contains 1000 realizations of convergence maps with the source redshift distribution of the CFHTLens survey; the simulated foregrounds are available for 91 different cosmological parameter variations of the triplet :math:`(\Omega_m,w,\sigma_8)`

This is an example on how you can use the LensTools API to interact with the simulations: suppose you have a local copy of the IGS1 simulations, which you saved in '/user/igs1' and you didn't modify the original directory tree. Then here is how you can interact with the maps

::

	>>> from lenstools.simulations import IGS1

	#Instantiate one of the models, for example the reference one
	>>> igs1_ref = IGS1(root_path="/user/igs1",H0=72.0,Om0=0.26,sigma8=0.798,ns=0.96,w0=-1)

	#Now you can access the names of the individual realizations
	>>> igs1_ref.getNames(2,z=1.0)
	
	/user/igs1/m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.798/Maps/WL-conv_m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.798_4096xy_0002r_0029p_0100z_og.gre.fit
	
	>>> igs1_ref.getNames(3,z=1.0,big_fiducial_set=True)
	
	/user/igs1/m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.798_f/Maps/WL-conv_m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.798_4096xy_0003r_0029p_0100z_og.gre.fit

	#Or you can load the images directly in memory
	>>> image = igs1_ref.load(2,z=1.0)

	#This wraps the image in a ConvergenceMap instance, now you can visualize the map
	>>> image.visualize()

.. figure:: ../../examples/conv_map.png

Read more in the :doc:`code` section to learn what you can do with ConvergenceMap instances; while LensTools provides the API to interact with a local copy of these simulated datasets, it does not provide the datasets themselves; to obtain the IGS1 dataset please email `Jan M. Kratochvil <jan.m.kratochvil@gmail.com>`_, to obtain the CFHTemu1 dataset please email `Andrea Petri <apetri@phys.columbia.edu>`_