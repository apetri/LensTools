import sys

####################################################
#######LensTools functionality######################
####################################################

from lenstools import ConvergenceMap,Ensemble,GaussianNoiseGenerator
from lenstools.index import PDF,Indexer
from lenstools.defaults import load_fits_default_convergence
from lenstools.simulations import IGS1

#####################################################

import numpy as np
import matplotlib.pyplot as plt

from emcee.utils import MPIPool

import logging

#########################################################################################
#########This function gets called on every map image and computes the histograms########
#########################################################################################

def compute_histograms(args):

	assert "map_id" in args.keys()
	assert "simulation_set" in args.keys()
	assert "smoothing_scales" in args.keys()
	assert "index" in args.keys()
	assert "generator" in args.keys()
	assert "bin_edges" in args.keys()

	assert len(args["index"].descriptor_list) == len(args["smoothing_scales"])

	z = 1.0

	#Get map name to analyze
	map_name = args["simulation_set"].getNames(z=z,realizations=[args["map_id"]])[0]

	#Load the convergence map
	convergence_map = ConvergenceMap.fromfilename(map_name,loader=load_fits_default_convergence)

	#Generate the shape noise map
	noise_map = args["generator"].getShapeNoise(z=z,ngal=15.0,seed=args["map_id"])

	#Add the noise
	convergence_map += noise_map

	#Measure the features
	hist_output = np.zeros(args["index"].size)
	for n,descriptor in enumerate(args["index"].descriptor_list):

		logging.debug("Processing {0} x {1} arcmin".format(map_name,args["smoothing_scales"][n]))

		smoothed_map = convergence_map.smooth(args["smoothing_scales"][n])
		v,hist_output[descriptor.first:descriptor.last] = smoothed_map.pdf(args["bin_edges"])

	#Return the histograms in array format
	return hist_output

############################################################
########################Main loop###########################
############################################################

if __name__=="__main__":
	
	#Initialize MPI pool
	try: 
		pool = MPIPool()
	except ValueError:
		pool = None
	
	#Set logging level
	logging.basicConfig(level=logging.DEBUG)
	
	if (pool is not None) and not(pool.is_master()):
	
		pool.wait()
		sys.exit(0)
	
	#Root path of IGS1 maps
	root_path = "/Users/andreapetri/Documents/Columbia/spurious_shear/convergence_maps"
	
	#Smoothing scales in arcmin
	smoothing_scales = [0.1,0.5,1.0,2.0]
	bin_edges = np.ogrid[-0.15:0.15:128j]
	bin_midpoints = 0.5*(bin_edges[1:] + bin_edges[:-1])
	
	#Create smoothing scale index for the histogram
	idx = Indexer.stack([PDF(bin_edges) for scale in smoothing_scales])
	
	#Create IGS1 simulation set object to look for the right simulations
	simulation_set = IGS1(root_path=root_path)
	
	#Look at a sample map
	sample_map = ConvergenceMap.fromfilename(simulation_set.getNames(z=1.0,realizations=[1])[0],loader=load_fits_default_convergence)
	
	#Initialize Gaussian shape noise generator
	generator = GaussianNoiseGenerator.forMap(sample_map)
	
	#Build Ensemble instance with the maps to analyze
	map_ensemble = Ensemble.fromfilelist(range(1,4))
	
	#Measure the histograms and load the data in the ensemble
	map_ensemble.load(callback_loader=compute_histograms,pool=pool,simulation_set=simulation_set,smoothing_scales=smoothing_scales,index=idx,generator=generator,bin_edges=bin_edges)
	
	if pool is not None:
		pool.close()

	##########################################################################################################################################
	###############################Ensemble data available at this point for covariance, PCA, etc...##########################################
	##########################################################################################################################################
	
	#Plot results to check
	fig,ax = plt.subplots(len(smoothing_scales),1)
	for i in range(len(smoothing_scales)):
		
		for n in range(map_ensemble.num_realizations):
			ax[i].plot(bin_midpoints,map_ensemble.data[n][idx[i].first:idx[i].last])
	
		ax[i].set_xlabel(r"$\kappa$")
		ax[i].set_ylabel(r"$P(\kappa)$")
		ax[i].set_title(r"${0:.1f}^\prime={1:.1f}$pix".format(smoothing_scales[i],smoothing_scales[i] * sample_map.kappa.shape[0]/(sample_map.side_angle*60.0)))
	
	
	fig.tight_layout()
	fig.savefig("histograms.png")

