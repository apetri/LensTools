from lenstools.statistics.ensemble import Ensemble
from lenstools.utils.defaults import default_callback_loader,peaks_loader
from lenstools.utils.decorators import Parallelize

import logging

import numpy as np
import matplotlib.pyplot as plt

logging.basicConfig(level=logging.DEBUG)

@Parallelize.masterworker
def main(pool): 

	l_edges = np.arange(200.0,50000.0,200.0)
	l = 0.5*(l_edges[:-1] + l_edges[1:])

	conv_ensemble = Ensemble.compute(["Data/conv1.fit","Data/conv2.fit","Data/conv3.fit","Data/conv4.fit"],callback_loader=default_callback_loader,pool=pool,l_edges=l_edges)

	fig,ax = plt.subplots()
	for n in range(len(conv_ensemble)):
		ax.plot(l,l*(l+1)*conv_ensemble.iloc[n]/(2.0*np.pi),label="Map {0}".format(n+1),linestyle="--")

	mean = conv_ensemble.mean(0)
	errors = np.sqrt(conv_ensemble.covariance().values.diagonal())

	ax.errorbar(l,l*(l+1)*mean/(2.0*np.pi),yerr=l*(l+1)*errors/(2.0*np.pi),label="Mean")

	ax.set_xscale("log")
	ax.set_yscale("log")
	ax.set_xlabel(r"$l$")
	ax.set_ylabel(r"$l(l+1)P_l/2\pi$")
	ax.legend(loc="upper left")

	plt.savefig("power_ensemble.png")

if __name__=="__main__":
	main(None)