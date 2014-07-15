import sys

try: 

	from lenstools import Ensemble
	from lenstools.index import Indexer,PowerSpectrum,Peaks
	from lenstools.defaults import convergence_measure_all

except ImportError:

	sys.path.append("..")
	from lenstools import Ensemble
	from lenstools.index import Indexer,PowerSpectrum,Peaks
	from lenstools.defaults import convergence_measure_all

import logging

import numpy as np
import matplotlib.pyplot as plt

logging.basicConfig(level=logging.DEBUG)

map_list = ["conv1.fit","conv2.fit","conv3.fit","conv4.fit"]
l_edges = np.arange(200.0,50000.0,200.0)
thresholds_pk = np.arange(-1.0,5.0,0.2)

def test_index():

	#Decide the statistical descriptors to measure, and build an index
	idx = Indexer.stack([PowerSpectrum(l_edges),Peaks(thresholds_pk)])
	l = idx[0].l
	v = idx[1].midpoints

	#Initiate the statistical ensemble
	ens = Ensemble.fromfilelist(map_list)

	#Load measurements into the ensemble
	ens.load(callback_loader=convergence_measure_all,pool=None,index=idx)

	#Split the ensemble in power_spectrum and peaks
	ens_pow,ens_pk = ens.split(idx)

	#####################################################################

	#Plot to check
	fig,ax = plt.subplots(1,2,figsize=(16,8))
	for i in range(ens.num_realizations):
		
		ax[0].plot(l,l*(l+1)*ens_pow.data[i]/(2.0*np.pi))
		ax[1].plot(v,ens_pk.data[i])

	ax[0].set_xscale("log")
	ax[0].set_yscale("log")

	ax[0].set_xlabel(r"$l$")
	ax[0].set_ylabel(r"$l(l+1)P_l/2\pi$")

	ax[1].set_xlabel(r"$\nu$")
	ax[1].set_ylabel(r"$dN/d\nu$")

	fig.tight_layout()

	plt.savefig("conv_all.png")
	plt.clf()


