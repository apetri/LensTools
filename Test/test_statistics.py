import sys,os

try:
	
	from lenstools import Ensemble
	from lenstools.defaults import default_callback_loader

except ImportError:

	sys.path.append("..")
	from lenstools import Ensemble
	from lenstools.defaults import default_callback_loader

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

map_list = ["conv1.fit","conv2.fit","conv3.fit","conv4.fit"]
l_edges = np.arange(200.0,50000.0,200.0)
l = 0.5*(l_edges[:-1] + l_edges[1:])

conv_ensemble = Ensemble.fromfilelist(map_list)
conv_ensemble.load(callback_loader=default_callback_loader,pool=pool,l_edges=l_edges)

if pool is not None:
	pool.close()

def test_shape():

	assert conv_ensemble.num_realizations==len(map_list)
	assert conv_ensemble.data.shape==(len(map_list),len(l_edges)-1)

def test_power_plot():

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
	plt.clf()

def test_chi2():

	conv_ensemble1 = Ensemble.fromfilelist(map_list[0:2])
	conv_ensemble1.load(callback_loader=default_callback_loader,pool=None,l_edges=l_edges)

	print("chi2 difference = {0}".format(conv_ensemble - conv_ensemble1))