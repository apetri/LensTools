import sys,os
	
from .. import Ensemble
from ..utils.defaults import measure_power_spectrum,peaks_loader

try:

	from emcee.utils import MPIPool
	MPIPool = MPIPool

except ImportError:

	MPIPool = None

import logging

from .. import dataExtern

import numpy as np
import pandas as pd
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
except (TypeError,ImportError):
	pool = None

#The only parallelized part is the loading of the ensemble (that's the computationally expensive part)

if (pool is not None) and not(pool.is_master()):

	pool.wait()
	sys.exit(0)

map_list = [os.path.join(dataExtern(),"conv1.fit"),os.path.join(dataExtern(),"conv2.fit"),os.path.join(dataExtern(),"conv3.fit"),os.path.join(dataExtern(),"conv4.fit")]
l_edges = np.arange(200.0,50000.0,200.0)
thresholds_pk = np.arange(-1.0,5.0,0.2)

l = 0.5*(l_edges[:-1] + l_edges[1:])

conv_ensemble = Ensemble.compute(file_list=map_list,callback_loader=measure_power_spectrum,pool=pool,l_edges=l_edges,columns=pd.Index(l,name="ell"))

if pool is not None:
	pool.close()

def test_shape():

	assert conv_ensemble.nobs==len(map_list)
	assert conv_ensemble.shape==(len(map_list),len(l_edges)-1)

def test_power_plot():

	fig,ax = plt.subplots()
	for n in range(len(conv_ensemble)):
		ax.plot(l,l*(l+1)*conv_ensemble.iloc[n].values/(2.0*np.pi),label="Map {0}".format(n+1),linestyle="--")

	mean = conv_ensemble.mean(0).values
	errors = np.sqrt(conv_ensemble.covariance().values.diagonal())

	ax.errorbar(l,l*(l+1)*mean/(2.0*np.pi),yerr=l*(l+1)*errors/(2.0*np.pi),label="Mean")

	ax.set_xscale("log")
	ax.set_yscale("log")
	ax.set_xlabel(r"$l$")
	ax.set_ylabel(r"$l(l+1)P_l/2\pi$")
	ax.legend(loc="upper left")

	plt.savefig("power_ensemble.png")
	plt.clf()

def test_chi2():
	conv_ensemble1 = Ensemble.compute(file_list=map_list[0:2],callback_loader=measure_power_spectrum,pool=None,l_edges=l_edges,columns=pd.Index(l,name="ell"))
	print("chi2 difference = {0}".format(conv_ensemble.compare(conv_ensemble1)))

def test_pca():

	pca_ensemble = Ensemble.read(os.path.join(dataExtern(),"ensemble_pca.npy"))
	pca = pca_ensemble.principalComponents()
	assert len(pca.explained_variance_)==pca_ensemble.shape[1]
	
	fig,ax = plt.subplots(1,2,figsize=(16,8))
	ax[0].plot(pca.explained_variance_)
	ax[1].plot(pca.explained_variance_.cumsum())
	ax[0].set_xlabel(r"$n$")
	ax[1].set_xlabel(r"$n$")
	ax[0].set_ylabel(r"$\lambda_n$")
	ax[1].set_ylabel(r"$\sum^n\lambda_n$")

	fig.savefig("pca.png")

def test_add():

	conv_ensemble1 = Ensemble.compute(file_list=map_list[0:2],callback_loader=measure_power_spectrum,pool=None,l_edges=l_edges,columns=pd.Index(l,name="ell"))
	conv_ensemble2 = Ensemble.compute(file_list=map_list[2:],callback_loader=measure_power_spectrum,pool=None,l_edges=l_edges,columns=pd.Index(l,name="ell"))
	conv_ensemble_union = Ensemble.concat([conv_ensemble1,conv_ensemble2],axis=0,ignore_index=True)

	assert conv_ensemble_union.nobs == 4
	assert conv_ensemble_union.shape[0] == 4
	assert conv_ensemble_union.shape[1] == conv_ensemble1.shape[1]

def test_multiply():

	conv_ensemble_peaks = Ensemble.compute(file_list=map_list,callback_loader=peaks_loader,pool=None,thresholds=thresholds_pk)
	conv_ensemble_both = Ensemble.concat([conv_ensemble,conv_ensemble_peaks],axis=1)

	assert conv_ensemble_both.nobs == 4
	assert conv_ensemble_both.shape[0] == 4
	assert conv_ensemble_both.shape[1] == len(l_edges) + len(thresholds_pk) - 2

def test_save_and_load():

	conv_ensemble.save("ensemble_saved.npy")
	conv_ensemble.save("ensemble_saved",format="matlab",appendmat=True)
	conv_ensemble_new = Ensemble.read("ensemble_saved.npy")

	assert conv_ensemble_new.nobs == conv_ensemble.nobs
	assert conv_ensemble_new.shape == conv_ensemble.shape

def test_group():

	conv_ensemble_sparse = Ensemble.compute(file_list=map_list,callback_loader=measure_power_spectrum,pool=pool,l_edges=l_edges)
	conv_ensemble_sparse = conv_ensemble_sparse.group(group_size=2,kind="sparse").mean()
	
	assert conv_ensemble_sparse.nobs==2

	conv_ensemble_contiguous = Ensemble.compute(file_list=map_list,callback_loader=measure_power_spectrum,pool=pool,l_edges=l_edges)
	conv_ensemble_contiguous = conv_ensemble_contiguous.group(group_size=2,kind="contiguous").mean()
	
	assert conv_ensemble_contiguous.nobs==2

	fig,ax = plt.subplots()
	for n in range(conv_ensemble.nobs):
		ax.plot(l,l*(l+1)*conv_ensemble.values[n]/(2.0*np.pi),label="Original {0}".format(n+1),linestyle="-")

	for n in range(conv_ensemble_sparse.nobs):
		ax.plot(l,l*(l+1)*conv_ensemble_sparse.values[n]/(2.0*np.pi),label="Sparse {0}".format(n+1),linestyle="--")

	for n in range(conv_ensemble_contiguous.nobs):
		ax.plot(l,l*(l+1)*conv_ensemble_contiguous.values[n]/(2.0*np.pi),label="Contiguous {0}".format(n+1),linestyle="-.")

	ax.set_xscale("log")
	ax.set_yscale("log")
	ax.set_xlabel(r"$l$")
	ax.set_ylabel(r"$l(l+1)P_l/2\pi$")
	ax.legend(loc="upper left",prop={"size":7})

	plt.savefig("power_ensemble_grouped.png")
	plt.clf()


def test_subset():

	conv_subset = conv_ensemble.iloc[range(2)]
	assert conv_subset.nobs==2

	fig,ax = plt.subplots()
	ax.plot(l,l*(l+1)*conv_subset.values[0]/(2.0*np.pi),label="1")
	ax.plot(l,l*(l+1)*conv_subset.values[1]/(2.0*np.pi),label="2")

	conv_subset = conv_ensemble.iloc[range(2,4)]
	assert conv_subset.nobs==2

	ax.plot(l,l*(l+1)*conv_subset.values[0]/(2.0*np.pi),label="3")
	ax.plot(l,l*(l+1)*conv_subset.values[1]/(2.0*np.pi),label="4")

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
	l_cut = list(filter(lambda ell:ell>=10000.0 and ell<=30000.0,conv_ensemble.columns))
	conv_ensemble_cut = conv_ensemble[l_cut]
	assert conv_ensemble_cut.shape[1] == len(l_cut)
	l_cut = np.array(l_cut)

	#Plot
	ax.plot(l_cut,l_cut*(l_cut+1)*conv_ensemble_cut.mean(0).values/(2.0*np.pi),label="Cut",color="yellow")

	ax.set_xscale("log")
	ax.set_yscale("log")
	ax.set_xlabel(r"$l$")
	ax.set_ylabel(r"$l(l+1)P_l/2\pi$")
	ax.legend(loc="upper left")

	plt.savefig("power_ensemble_cut.png")
	plt.clf()


def test_differentiate():

	conv_ensemble_peaks_cumulative = Ensemble.compute(file_list=map_list,callback_loader=peaks_loader,pool=None,thresholds=thresholds_pk).cumsum(1)
	diff_ensemble = conv_ensemble_peaks_cumulative.apply(lambda s:s.diff(),axis=1)

	fig,ax = plt.subplots()
	for i in range(diff_ensemble.nobs):
		ax.plot(0.5*(thresholds_pk[:-1]+thresholds_pk[1:]),diff_ensemble.values[i])
		
	ax.set_xlabel(r"$\kappa$")
	ax.set_ylabel(r"$P(\kappa)$")

	fig.savefig("ensemble_differentiate.png")


def test_selfChi2():

	ens = Ensemble.read(os.path.join(dataExtern(),"all","Om0.295_Ol0.705_w-1.878_ns0.960_si0.100","subfield1","sigma05","power_spectrum.npy"))
	chi2 = ens.selfChi2()
	assert chi2.shape[0]==ens.shape[0]

	#Plot histogram
	fig,ax = plt.subplots()
	n,bins,patch = ax.hist(chi2.values,bins=50,normed=True,histtype="stepfilled",alpha=0.5)

	#Compare to chi2 distribution
	ax.plot(stats.chi2.pdf(bins,ens.shape[1]))

	#Labels
	ax.set_xlabel(r"$\chi^2$")
	ax.set_ylabel(r"$P(\chi^2)$")

	#Save figure
	fig.savefig("self_chi2.png")



