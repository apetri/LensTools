try:
	
	from lenstools import Ensemble
	from lenstools.defaults import default_callback_loader

except ImportError:
	
	import sys
	sys.path.append("..")
	from lenstools import Ensemble
	from lenstools.defaults import default_callback_loader

import numpy as np
import matplotlib.pyplot as plt

map_list = ["conv1.fit","conv2.fit","conv3.fit","conv4.fit"]
l_edges = np.arange(200.0,50000.0,200.0)
l = 0.5*(l_edges[:-1] + l_edges[1:])

conv_ensemble = Ensemble.fromfilelist(map_list,callback_loader=default_callback_loader,l_edges=l_edges)

def test_shape():

	assert conv_ensemble.num_realizations==len(map_list)
	assert conv_ensemble.data.shape==(len(map_list),len(l_edges)-1)

def test_power_plot():

	fig,ax = plt.subplots()
	for n in range(conv_ensemble.num_realizations):
		ax.plot(l,l*(l+1)*conv_ensemble.data[n]/(2.0*np.pi),label="Map {0}".format(n+1))

	ax.set_xscale("log")
	ax.set_yscale("log")
	ax.set_xlabel(r"$l$")
	ax.set_ylabel(r"$l(l+1)P_l/2\pi$")
	ax.legend(loc="upper left")

	plt.savefig("power_ensemble.png")
	plt.clf()