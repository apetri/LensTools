"""A collection of tools widely used in Weak Gravitational Lensing data analyses

.. moduleauthor:: Andrea Petri <apetri@phys.columbia.edu>


"""

__version__ = "1.2"

import sys

from .utils.configuration import configuration

from .catalog import Catalog,ShearCatalog,FlexionCatalog,MockCatalog 
from .image.convergence import ConvergenceMap,OmegaMap,Mask,CMBTemperatureMap
from .image.shear import ShearMap
from .image.flexion import FlexionMap
from .image.noise import GaussianNoiseGenerator
from .statistics.ensemble import Ensemble

from .statistics import contours,constraints

from .pipeline.simulation import SimulationBatch

from .legacy import index

import os,pkg_resources

if sys.version_info.major>=3:
	from urllib import request
else:
	import urllib2 as request

import tarfile

#External data needed by tests, may be downloaded
data_directory = os.getenv("LENSTOOLS_DATA")
if data_directory is None:
	data_directory = "LT_Data"

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
	
	data_url = "https://www.dropbox.com/s/oh526rkrhiy3m8u/data.tar.gz"
	data_filename = os.path.join(path,"data.tar.gz")

	#Download the file
	response = request.urlopen(data_url)
	with open(data_filename,"wb") as datafile:
		print("[+] Downloading {0}...".format(data_url))
		datafile.write(response.read())

	#Unpack the archive
	print("[+] Unpacking data archive {0}".format(data_filename))
	with tarfile.open(data_filename,"r:gz") as tar:
def is_within_directory(directory, target):
	
	abs_directory = os.path.abspath(directory)
	abs_target = os.path.abspath(target)

	prefix = os.path.commonprefix([abs_directory, abs_target])
	
	return prefix == abs_directory

def safe_extract(tar, path=".", members=None, *, numeric_owner=False):

	for member in tar.getmembers():
		member_path = os.path.join(path, member.name)
		if not is_within_directory(path, member_path):
			raise Exception("Attempted Path Traversal in Tar File")

	tar.extractall(path, members, numeric_owner=numeric_owner) 
	

safe_extract(tar, path)

	#Remove the archive
	os.remove(data_filename)

	#Make the necessary rename
	os.rename(os.path.join(path,"Data"),os.path.join(path,os.path.basename(data_directory)))
	print("[*] Data directory is available at {0}".format(dataExtern()))




