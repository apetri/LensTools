import sys

try: 

	from lenstools import Ensemble
	from lenstools.index import Indexer,PowerSpectrum,PDF,Peaks,MinkowskiAll,Moments
	from lenstools.defaults import convergence_measure_all

except ImportError:

	sys.path.append("..")
	from lenstools import Ensemble
	from lenstools.index import Indexer,PowerSpectrum,PDF,Peaks,MinkowskiAll,Moments
	from lenstools.defaults import convergence_measure_all

import logging

import numpy as np
import matplotlib.pyplot as plt

logging.basicConfig(level=logging.DEBUG)

map_list = ["Data/conv1.fit","Data/conv2.fit","Data/conv3.fit","Data/conv4.fit"]
l_edges = np.arange(200.0,50000.0,200.0)
thresholds_pk = np.arange(-1.0,5.0,0.2)
thresholds_mf = np.arange(-2.0,2.0,0.2)

def test_index():

	#Decide the statistical descriptors to measure, and build an index
	idx = Indexer.stack([PowerSpectrum(l_edges),Peaks(thresholds_pk,norm=True),MinkowskiAll(thresholds_mf,norm=True),PDF(thresholds_mf,norm=True),Moments(connected=True)])
	l = idx[0].l
	v = idx[1].midpoints
	v_mf = idx[2].midpoints

	#Initiate the statistical ensemble
	ens = Ensemble.fromfilelist(map_list)

	#Load measurements into the ensemble (this is the expensive part!!!)
	ens.load(callback_loader=convergence_measure_all,pool=None,index=idx)

	#Split the ensemble in power_spectrum,peaks, and the second and third minkowski functional
	mink_idx = idx[2].separate()
	subset_idx = Indexer([idx[0],idx[1],idx[3],mink_idx[2],idx[-1]])

	ens_pow,ens_pk,ens_pdf,ens_mink2,ens_mom = ens.split(subset_idx)

	#####################################################################

	#Plot to check
	fig,ax = plt.subplots(2,2,figsize=(16,16))
	for i in range(ens.num_realizations):
		
		ax[0,0].plot(l,l*(l+1)*ens_pow.data[i]/(2.0*np.pi))
		ax[0,1].plot(v,ens_pk.data[i])
		ax[1,0].plot(v_mf,ens_pdf.data[i])
		ax[1,1].plot(v_mf,ens_mink2.data[i])

	ax[0,0].set_xscale("log")
	ax[0,0].set_yscale("log")

	ax[0,0].set_xlabel(r"$l$")
	ax[0,0].set_ylabel(r"$l(l+1)P_l/2\pi$")

	ax[0,1].set_xlabel(r"$\nu$")
	ax[0,1].set_ylabel(r"$dN/d\nu$")

	ax[1,0].set_xlabel(r"$\nu$")
	ax[1,0].set_ylabel(r"$P(\nu)$")

	ax[1,1].set_xlabel(r"$\nu$")
	ax[1,1].set_ylabel(r"$V_2(\nu)$")

	fig.tight_layout()

	plt.savefig("conv_all.png")
	plt.clf()

	#Save moments to check
	np.savetxt("moments.txt",ens_mom.mean())


