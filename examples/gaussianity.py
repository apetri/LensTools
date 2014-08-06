"""
Tests the gaussianity of simulated noise maps by measuring their cubic and quartic moments

"""
import sys

from lenstools import ConvergenceMap,Ensemble,GaussianNoiseGenerator
from lenstools.defaults import load_fits_default_convergence
from lenstools.simulations import IGS1

import numpy as np

import logging

from emcee.utils import MPIPool


def generate_and_measure(map_id,generator,power_func):

	#Generate the noise map
	logging.debug("Processing mock map {0}".format(map_id))
	conv_map = generator.fromConvPower(power_func=power_func,seed=map_id,bounds_error=False,fill_value=0.0)

	#Measure its moments
	return conv_map.moments(connected=True,dimensionless=True)

def measure_from_IGS1(filename):

	#Read in the map
	logging.debug("Processing IGS1 map {0}".format(filename))
	conv_map = ConvergenceMap.fromfilename(filename,loader=load_fits_default_convergence)

	#Smooth 1 arcmin
	conv_map.smooth(1.0,inplace=True)

	#Measure the moments
	return conv_map.moments(connected=True,dimensionless=True)




logging.basicConfig(level=logging.DEBUG)

try: 
	pool = MPIPool()
except ValueError:
	pool = None

if (pool is not None) and not(pool.is_master()):

	pool.wait()
	sys.exit(0)



map_mock_ids = range(int(sys.argv[1]))

igs1_set = IGS1(root_path="/astro/astronfs01/jank/Storage/wl/IG/m-series")
map_igs1_ids = igs1_set.getNames(z=2.0,realizations=range(1,int(sys.argv[1])+1))

gen = GaussianNoiseGenerator(shape=(2048,2048),side_angle=3.41,label="convergence") 
power_func = np.loadtxt("Data/ee4e-7.txt",unpack=True)

ens_mock = Ensemble.fromfilelist(map_mock_ids)
ens_igs1 = Ensemble.fromfilelist(map_igs1_ids)

ens_mock.load(callback_loader=generate_and_measure,pool=pool,generator=gen,power_func=power_func)
ens_igs1.load(callback_loader=measure_from_IGS1,pool=pool)

if pool is not None:
	pool.close()

np.savetxt("moments_mock.txt",np.array([ens_mock.mean(),np.sqrt(ens_mock.covariance().diagonal())]))
np.savetxt("moments_igs1.txt",np.array([ens_igs1.mean(),np.sqrt(ens_igs1.covariance().diagonal())]))

logging.info("Done!")
