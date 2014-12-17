import sys

try:
	
	from lenstools import Ensemble
	from lenstools.defaults import default_callback_loader,peaks_loader,convergence_measure_all
	from lenstools.index import Indexer,MinkowskiAll

except ImportError:

	sys.path.append("..")
	from lenstools import Ensemble
	from lenstools.defaults import default_callback_loader,peaks_loader,convergence_measure_all
	from lenstools.index import Indexer,MinkowskiAll

try:

	from emcee.utils import MPIPool
	MPIPool = MPIPool

except ImportError:

	MPIPool = None

import logging

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

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
		ax.plot(l,l*(l+1)*conv_ensemble[n]/(2.0*np.pi),label="Map {0}".format(n+1),linestyle="--")

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

def test_pca():

	pca_ensemble = Ensemble.read("Data/ensemble_pca.npy")
	pca = pca_ensemble.principalComponents()
	assert len(pca.explained_variance_)==pca_ensemble.data.shape[1]
	
	fig,ax = plt.subplots(1,2,figsize=(16,8))
	ax[0].plot(pca.explained_variance_)
	ax[1].plot(pca.explained_variance_.cumsum())
	ax[0].set_xlabel(r"$n$")
	ax[1].set_xlabel(r"$n$")
	ax[0].set_ylabel(r"$\lambda_n$")
	ax[1].set_ylabel(r"$\sum^n\lambda_n$")

	fig.savefig("pca.png")

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
	conv_ensemble.save("ensemble_saved",format="matlab",appendmat=True)
	conv_ensemble_new = Ensemble.read("ensemble_saved.npy")

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

def test_subset():

	conv_subset = conv_ensemble.subset(range(2))
	assert conv_subset.num_realizations==2

	fig,ax = plt.subplots()
	ax.plot(l,l*(l+1)*conv_subset[0]/(2.0*np.pi),label="1")
	ax.plot(l,l*(l+1)*conv_subset[1]/(2.0*np.pi),label="2")

	conv_subset = conv_ensemble.subset(range(2,4))
	assert conv_subset.num_realizations==2

	ax.plot(l,l*(l+1)*conv_subset[0]/(2.0*np.pi),label="3")
	ax.plot(l,l*(l+1)*conv_subset[1]/(2.0*np.pi),label="4")

	ax.set_xscale("log")
	ax.set_yscale("log")
	ax.set_xlabel(r"$l$")
	ax.set_ylabel(r"$l(l+1)P_l/2\pi$")
	ax.legend(loc="upper left")

	fig.savefig("power_ensemble_subset.png")

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


def test_differentiate():

	thresholds = np.arange(-0.04,0.12,0.001)
	midpoints = 0.5*(thresholds[:-1] + thresholds[1:])

	index = Indexer.stack([MinkowskiAll(thresholds)])
	index_separate = Indexer(MinkowskiAll(thresholds).separate())
	
	diff_ensemble = Ensemble.fromfilelist(map_list)
	diff_ensemble.load(callback_loader=convergence_measure_all,index=index)

	ensemble_0 = diff_ensemble.split(index_separate)[0]
	ensemble_pdf = ensemble_0.differentiate(step=thresholds[0]-thresholds[1])

	fig,ax = plt.subplots()
	for i in range(ensemble_0.num_realizations):
		ax.plot(0.5*(midpoints[:-1]+midpoints[1:]),ensemble_pdf[i])
		
	ax.set_xlabel(r"$\kappa$")
	ax.set_ylabel(r"$P(\kappa)$")

	fig.savefig("ensemble_differentiate.png")


def test_selfChi2():

	ens = Ensemble.read("Data/all/Om0.295_Ol0.705_w-1.878_ns0.960_si0.100/subfield1/sigma05/power_spectrum.npy")
	chi2 = ens.selfChi2()
	assert chi2.shape[0]==ens.data.shape[0]

	#Plot histogram
	fig,ax = plt.subplots()
	n,bins,patch = ax.hist(chi2,bins=50,normed=True,histtype="stepfilled",alpha=0.5)

	#Compare to chi2 distribution
	ax.plot(stats.chi2.pdf(bins,ens.data.shape[1]))

	#Labels
	ax.set_xlabel(r"$\chi^2$")
	ax.set_ylabel(r"$P(\chi^2)$")

	#Save figure
	fig.savefig("self_chi2.png")



