"""
Tests the gaussianity of simulated noise maps by measuring their cubic and quartic moments

"""

from lenstools import Ensemble,GaussianNoiseGenerator
import numpy as np

import logging

def generate_and_measure(args):

	assert "map_id" in args.keys()
	assert "generator" in args.keys()
	assert "power_func" in args.keys()

	#Generate the noise map
	logging.debug("Processing map {0}".format(args["map_id"]))
	conv_map = args["generator"].fromConvPower(power_func=args["power_func"],seed=args["map_id"],bounds_error=False,fill_value=0.0)

	#Measure its moments
	return conv_map.moments(connected=True,dimensionless=True)

logging.basicConfig(level=logging.DEBUG)

map_ids = range(4)

gen = GaussianNoiseGenerator(shape=(2048,2048),side_angle=3.41,label="convergence") 
power_func = np.loadtxt("ee4e-7.txt",unpack=True)

ens = Ensemble.fromfilelist(map_ids)

ens.load(callback_loader=generate_and_measure,pool=None,generator=gen,power_func=power_func)

np.savetxt("moments.txt",np.array([ens.mean(),np.sqrt(ens.covariance().diagonal())]))