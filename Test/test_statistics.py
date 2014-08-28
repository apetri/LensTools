import sys

try:
	
	from lenstools import Ensemble
	from lenstools.defaults import default_callback_loader,peaks_loader

except ImportError:

	sys.path.append("..")
	from lenstools import Ensemble
	from lenstools.defaults import default_callback_loader,peaks_loader

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
thresholds_pk = np.arange(-1.0,5.0,0.2)

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

	print("chi2 difference = {0}".format(conv_ensemble.compare(conv_ensemble1)))

def test_add():

	conv_ensemble1 = Ensemble.fromfilelist(map_list[0:2])
	conv_ensemble2 = Ensemble.fromfilelist(map_list[2:])

	conv_ensemble1.load(callback_loader=default_callback_loader,pool=None,l_edges=l_edges)
	conv_ensemble2.load(callback_loader=default_callback_loader,pool=None,l_edges=l_edges)

	conv_ensemble_union = conv_ensemble1 + conv_ensemble2

	assert conv_ensemble_union.num_realizations == 4
	assert len(conv_ensemble_union.file_list) == 4
	assert conv_ensemble_union.data.shape[0] == 4
	assert conv_ensemble_union.data.shape[1] == conv_ensemble1.data.shape[1]

def test_multiply():

	conv_ensemble_peaks = Ensemble.fromfilelist(map_list)
	conv_ensemble_peaks.load(callback_loader=peaks_loader,pool=None,thresholds=thresholds_pk)

	conv_ensemble_both = conv_ensemble * conv_ensemble_peaks

	assert conv_ensemble_both.num_realizations == 4
	assert conv_ensemble_both.data.shape[0] == 4
	assert conv_ensemble_both.data.shape[1] == len(l_edges) + len(thresholds_pk) - 2

def test_save_and_load():

	conv_ensemble.save("ensemble_saved.npy")
	conv_ensemble.savemat("ensemble_saved",appendmat=True)
	conv_ensemble_new = Ensemble.fromfilelist(["ensemble_saved.npy"])

	conv_ensemble_new.load(from_old=True)

	assert conv_ensemble_new.num_realizations == conv_ensemble.num_realizations
	assert conv_ensemble_new.data.shape == conv_ensemble.data.shape

def test_group():

	conv_ensemble_sparse = Ensemble.fromfilelist(map_list)
	conv_ensemble_sparse.load(callback_loader=default_callback_loader,pool=pool,l_edges=l_edges)
	conv_ensemble_sparse.group(group_size=2,kind="sparse")
	
	assert conv_ensemble_sparse.num_realizations==2

	conv_ensemble_contiguous = Ensemble.fromfilelist(map_list)
	conv_ensemble_contiguous.load(callback_loader=default_callback_loader,pool=pool,l_edges=l_edges)
	conv_ensemble_contiguous.group(group_size=2,kind="contiguous")
	
	assert conv_ensemble_contiguous.num_realizations==2

	fig,ax = plt.subplots()
	for n in range(conv_ensemble.num_realizations):
		ax.plot(l,l*(l+1)*conv_ensemble.data[n]/(2.0*np.pi),label="Original {0}".format(n+1),linestyle="-")

	for n in range(conv_ensemble_sparse.num_realizations):
		ax.plot(l,l*(l+1)*conv_ensemble_sparse.data[n]/(2.0*np.pi),label="Sparse {0}".format(n+1),linestyle="--")

	for n in range(conv_ensemble_contiguous.num_realizations):
		ax.plot(l,l*(l+1)*conv_ensemble_contiguous.data[n]/(2.0*np.pi),label="Contiguous {0}".format(n+1),linestyle="-.")

	ax.set_xscale("log")
	ax.set_yscale("log")
	ax.set_xlabel(r"$l$")
	ax.set_ylabel(r"$l(l+1)P_l/2\pi$")
	ax.legend(loc="upper left",prop={"size":7})

	plt.savefig("power_ensemble_grouped.png")
	plt.clf()

	return conv_ensemble_sparse._scheme,conv_ensemble_contiguous._scheme

def test_cut():

	fig,ax = plt.subplots()
	ax.plot(l,l*(l+1)*conv_ensemble.mean()/(2.0*np.pi),label="Full")

	#Perform the cut
	l_cut = conv_ensemble.cut(10000.0,30000.0,feature_label=l)
	assert conv_ensemble.data.shape[1] == len(l_cut)

	#Plot
	ax.plot(l_cut,l_cut*(l_cut+1)*conv_ensemble.mean()/(2.0*np.pi),label="Cut",color="yellow")

	ax.set_xscale("log")
	ax.set_yscale("log")
	ax.set_xlabel(r"$l$")
	ax.set_ylabel(r"$l(l+1)P_l/2\pi$")
	ax.legend(loc="upper left")

	plt.savefig("power_ensemble_cut.png")
	plt.clf()

