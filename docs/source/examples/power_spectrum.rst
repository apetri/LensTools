.. _power_spectrum::

Measure the power spectrum of a convergence map
===============================================

Snipped from lenstools.utils.defaults 

::

	################################################################################
	##########Default callback loader, loads in the measured power spectrum#########
	################################################################################

	def default_callback_loader(filename,l_edges):
		"""
	
		Default ensemble loader: reads a FITS data file containing a convergence map and measures its power spectrum

		:param args: A dictionary that contains all the relevant parameters as keys. Must have a "map_id" key
		:type args: Dictionary

		:returns: ndarray of the measured statistics

		:raises: AssertionError if the input dictionary doesn't have the required keywords

		"""

		logging.debug("Processing {0} power".format(filename))

		conv_map = ConvergenceMap.load(filename)
		l,Pl = conv_map.powerSpectrum(l_edges)
		return Pl

This is the actual code 

::

	from lenstools import Ensemble
	from lenstools.utils.defaults import default_callback_loader

	try:

		from emcee.utils import MPIPool
		MPIPool = MPIPool

	except ImportError:

		MPIPool = None

	import logging

	import numpy as np
	import matplotlib.pyplot as plt

	logging.basicConfig(level=logging.DEBUG)

	if MPIPool is None:
		logging.warning("You need to install emcee in order to test the parallel statistics features!!")

	try:
		logging.debug("Attempting to create MPIPool")
		pool = MPIPool()
		logging.debug("Succesfully created MPIPool!")
	except ValueError:
		logging.debug("No reason to create one, one process only!!")
		pool = None
	except TypeError:
		pool = None

	#The only parallelized part is the loading of the ensemble (that's the computationally expensive part)

	if (pool is not None) and not(pool.is_master()):

		pool.wait()
		sys.exit(0)

	map_list = ["Data/conv1.fit","Data/conv2.fit","Data/conv3.fit","Data/conv4.fit"]

	l_edges = np.arange(200.0,50000.0,200.0)
	l = 0.5*(l_edges[:-1] + l_edges[1:])

	conv_ensemble = Ensemble.fromfilelist(map_list)
	conv_ensemble.load(callback_loader=default_callback_loader,pool=pool,l_edges=l_edges)

	if pool is not None:
		pool.close()

	fig,ax = plt.subplots()
	for n in range(conv_ensemble.num_realizations):
		ax.plot(l,l*(l+1)*conv_ensemble.data[n]/(2.0*np.pi),label="Map {0}".format(n+1),linestyle="--")

	mean = conv_ensemble.mean()
	errors = np.sqrt(conv_ensemble.covariance().diagonal())

	ax.errorbar(l,l*(l+1)*mean/(2.0*np.pi),yerr=l*(l+1)*errors/(2.0*np.pi),label="Mean")

	ax.set_xscale("log")
	ax.set_yscale("log")
	ax.set_xlabel(r"$l$")
	ax.set_ylabel(r"$l(l+1)P_l/2\pi$")
	ax.legend(loc="upper left")

	plt.savefig("power_ensemble.png")

And this is the result 

.. figure:: ../../../examples/power_ensemble.png