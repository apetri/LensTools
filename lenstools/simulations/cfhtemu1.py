from __future__ import division,print_function,with_statement

import os,re
from pkg_resources import resource_filename

import numpy as np
from astropy.io import fits
from astropy.units import deg

from .igs1 import IGS1

######################################
#######Loader function for CFHT#######
######################################

def cfht_load(self,filename):

	kappa_file = fits.open(filename)
	angle = 3.4641016151377544 * deg
	kappa = kappa_file[0].data.astype(np.float)
	kappa_file.close()

	return angle,kappa

######################################
#######EMU1 simulations###############
######################################

class CFHTemu1(IGS1):

	"""
	Class handler of the weak lensing CFHTemu1 simulations set, inherits from IGS1; this simulation suite contains 91 different cosmological models based on 1 N-body simulation each. Each model has 1000 realizations for each of the 13 CFHT subfields 

	"""

	#Don't touch these! 
	_class_name = "CFHTemu1"
	_series_name = "emu1"
	_num_particles = 512
	_box_size_mpc = 240
	_data_loader = cfht_load

	@classmethod
	def getModels(cls,root_path="/default"):
		"""
		This class method uses a dictionary file to read in all the cosmological model parameters and instantiate the corresponding CFHTemu1 objects for each one of them

		:param root_path: path of your CFHT emu1 simulations copy
		:type root_path: str.

		:returns: list of CFHTemu1 instances

		"""

		#Build the complete filename
		dict_filename = resource_filename("lenstools",os.path.join("data","CFHTemu1.txt"))

		#Read in the dictionary file
		with open(dict_filename,"r") as dictfile:
			cosmo_id_strings = dictfile.readlines()

		#Read each cosmo identifier, squeeze out the cosmological parameters for each of them and instantiate the corresponding CFHTemu1 object
		model_list = list()
		
		for cosmo_id in cosmo_id_strings:

			#Squeeze out the cosmological parameters with a regular expression match
			m = re.match(r"Om([0-9.]{5})_Ol([0-9.]{5})_w([0-9\-.]+)_ns([0-9.]{5})_si([0-9.]{5})",cosmo_id)
			Om0,Ol0,w0,ns,sigma8 = m.groups()
			
			#Instantiate the model object
			model_list.append(cls(root_path=root_path,name=cosmo_id.rstrip("\n"),H0=70.0,Om0=float(Om0),w0=float(w0),sigma8=float(sigma8),ns=float(ns)))

		#Return the model list
		return model_list


	def getNames(self,realizations,subfield=1,smoothing=0.5):

		"""
		Get the full name of the CFHT emu1 maps, once a subfield and smoothing scale are specified

		:param subfield: the specific CFHT subfield you want to retrieve, must be between 1 and 13
		:type subfield: int.

		:param smoothing: smoothing scale of the maps you wish to retrieve, in arcmin
		:type smoothing: float.

		:param realizations: list of realizations to get the names of, the elements must be in [1,1000]
		:type realizations: list. or int.

		"""

		assert 1 <= subfield <= 13
		assert type(realizations) == list or type(realizations) == int

		#Build the file name
		root_path = self.root_path

		name = os.path.join(root_path,self._series_name + "-")
		name += self._box_string + "_" 
		name += self._cosmo_id_string 
		name = os.path.join(name,"subfield{0}".format(subfield)) 
		name = os.path.join(name,"sigma{0:02d}".format(int(smoothing*10)))
		name = os.path.join(name,"SIM_KS_sigma{0:02d}_subfield{1}_{2}-{3}_{4}_".format(int(smoothing*10),subfield,self._series_name,self._box_string,self._cosmo_id_string))

		#return the results
		if type(realizations) == int:
			return name + "{0:0004d}r.fit".format(realizations)
		else:
			return [name + "{0:0004d}r.fit".format(r) for r in realizations]


#########################################
########cfhtcov simulations##############
#########################################

class CFHTcov(CFHTemu1):

	"""
	Class handler of the weak lensing CFHTcov simulations set, inherits from CFHTemu1; this simulation suite contains 1000 realizations for each of the 13 CFHT subfields, based on 50 independent N-body simulations of a fiducial LambdaCDM universe. Useful to measure accurately descriptor covariance matrices

	"""

	#Don't touch these! 
	_class_name = "CFHTcov"
	_series_name = "cfhtcov"
	_num_particles = 512
	_box_size_mpc = 240


	@classmethod
	def getModels(cls,root_path="/default"):

		"""
		On call, this class method returns a CFHTcov instance initialized with the cosmological parameters of the only available model in the suite 

		:param root_path: path of your CFHTcov simulations copy
		:type root_path: str.

		:returns: CFHTcov instance initialized with the fiducial cosmological parameters

		"""

		return cls(root_path=root_path,name="Om{0:.3f}_Ol{1:.3f}_w{2:.3f}_ns{3:.3f}_si{4:.3f}".format(0.26,0.74,-1.0,0.960,0.800),H0=70.0,Om0=0.26,w0=-1.0,sigma8=0.800,ns=0.960)


	def getNames(self,realizations,subfield=1,smoothing=0.5):

		"""
		Get the full name of the CFHTcov maps, once a subfield and smoothing scale are specified

		:param subfield: the specific CFHT subfield you want to retrieve, must be between 1 and 13
		:type subfield: int.

		:param smoothing: smoothing scale of the maps you wish to retrieve, in arcmin
		:type smoothing: float.

		:param realizations: list of realizations to get the names of, the elements must be in [1,1000]
		:type realizations: list. or int.

		"""

		assert 1 <= subfield <= 13
		assert type(realizations) == list or type(realizations) == int

		#Build the file name
		root_path = self.root_path

		name = os.path.join(root_path,self._series_name + "-")
		name += self._box_string + "_" 
		name += self._cosmo_id_string 
		name = os.path.join(name,"subfield{0}".format(subfield)) 
		name = os.path.join(name,"sigma{0:02d}".format(int(smoothing*10)))
		name = os.path.join(name,"SIM_KS_sigma{0:02d}_subfield{1}_WL-only_{2}-{3}_{4}_{5}xy_".format(int(smoothing*10),subfield,self._series_name,self._box_string,self._cosmo_id_string,self._lens_plane_size))

		#return the results
		if type(realizations) == int:
			return name + "{0:0004d}r.fit".format(realizations)
		else:
			return [name + "{0:0004d}r.fit".format(r) for r in realizations]

