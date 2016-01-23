"""A collection of tools widely used in Weak Gravitational Lensing data analyses

.. moduleauthor:: Andrea Petri <apetri@phys.columbia.edu>


"""

__version__ = "0.5"

from lenstools.utils.configuration import configuration

from image.convergence import ConvergenceMap,Mask
from image.shear import ShearMap
from image.noise import GaussianNoiseGenerator
from .statistics.ensemble import Ensemble

import statistics.contours as contours
import statistics.constraints as constraints

from lenstools.pipeline.simulation import SimulationBatch

import legacy.index as index

import os,pkg_resources
import urllib2
import tarfile

#External data needed by tests, may be downloaded
data_directory = os.getenv("LENSTOOLS_DATA")
if data_directory is None:
	data_directory = "Data"

#Path to the data folder
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
	
	if not(os.path.exists(data_directory)):
		getTestData(os.path.dirname(os.path.abspath(data_directory)))

	return data_directory


#Download the test data into path
def getTestData(path="."):
	
	data_url = "http://www.plus-six.com/andrea/lenstools/data.tar.gz"
	data_filename = os.path.join(path,"data.tar.gz")

	#Download the file
	response = urllib2.urlopen(data_url)
	with open(data_filename,"wb") as datafile:
		print("[+] Downloading {0}...".format(data_url))
		datafile.write(response.read())

	#Unpack the archive
	with tarfile.open(data_filename,"r:gz") as tar:
		tar.extractall(path)

	#Remove the archive
	os.remove(data_filename)

	#Make the necessary rename
	os.rename(os.path.join(path,"Data"),os.path.join(path,os.path.basename(data_directory)))



