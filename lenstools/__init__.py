"""A collection of tools widely used in Weak Gravitational Lensing data analyses

.. moduleauthor:: Andrea Petri <apetri@phys.columbia.edu>


"""

__version__ = "0.5"

from image.convergence import ConvergenceMap,Mask
from image.shear import ShearMap
from image.noise import GaussianNoiseGenerator
from .statistics.ensemble import Ensemble

import statistics.index as index
import statistics.contours as contours
import statistics.constraints as constraints


#Path to the data folder
import os,pkg_resources

def data(name=None):

	if name is not None:
		
		full_path = pkg_resources.resource_filename("lenstools",os.path.join("data",name))
		if os.path.isfile(full_path):
			return full_path
		else:
			raise IOError("The file {0} does not exist!".format(full_path))

	else:

		#If no name provided just list all available resources
		full_path = pkg_resources.resource_filename("lenstools","data")
		return os.listdir(full_path)


def showData(name):
	
	path = data(name)
	with open(path,"r") as datafile:
		print(datafile.read())

def dataExtern():
	return "Data"