from __future__ import division,print_function,with_statement

import os,re

import numpy as np

from astropy.cosmology import FlatwCDM
from astropy.io import fits
from astropy.units import deg

from ..image.convergence import ConvergenceMap


######################################
#######Loader function for IGS1#######
######################################

def igs1_load(self,filename):

	kappa_file = fits.open(filename)
	angle = kappa_file[0].header["ANGLE"] * deg
	kappa = kappa_file[0].data.astype(np.float)
	kappa_file.close()

	return angle,kappa


######################################
###########IGS1 class#################
######################################

class IGS1(FlatwCDM):

	"""
	Class handler of the IGS1 simulations set, inherits the cosmological parameters from the astropy.cosmology.FlatwCDM class; the default parameter values are the fiducial ones

	"""

	#Don't touch these! 
	_class_name = "IGS1"
	_series_name = "m"
	_num_particles = 512
	_box_size_mpc = 240
	_lens_plane_size = 4096
	_data_loader = igs1_load


	def __init__(self,H0=72.0,Om0=0.26,w0=-1.0,sigma8=0.798,ns=0.960,root_path=None,name=None):

		super(IGS1,self).__init__(H0,Om0,w0=w0,name=name)
		self.sigma8 = sigma8
		self.ns = ns

		assert root_path is not None,"You must specify the root path of your {0} local copy!".format(self._class_name)

		self.root_path = root_path

		#Don't touch these! 
		self._cosmo_id_string =  "Om{0:.3f}_Ol{1:.3f}_w{2:.3f}_ns{3:.3f}_si{4:.3f}".format(self.Om0,1.0-self.Om0,self.w0,self.ns,self.sigma8)
		self._box_string = str(self._num_particles)+"b"+str(self._box_size_mpc)
		self._full_path = os.path.join(self.root_path,self._series_name+"-"+self._box_string+"_"+self._cosmo_id_string)

	def __repr__(self):

		astropy_string = super(IGS1,self).__repr__()
		pieces = astropy_string.split(",")
		si8_piece = u" sigma8={0}".format(self.sigma8)
		ns_piece = u" ns={0}".format(self.ns)

		return ",".join(pieces[:3] + [si8_piece,ns_piece] + pieces[3:])

	def _plane_id(self,z):

		if z==1.0:
			return "0029p"
		elif z==1.5:
			return "0038p"
		elif z==2.0:
			return "0046p"
		else:
			raise ValueError("IGS1 doesn't have maps at redshift {0}".format(z))

	def squeeze(self,with_ns=False):
		
		"""
		Returns the cosmological parameters of the model in numpy array form

		:param with_ns: if True returns also ns as the last parameter
		:type with_ns: bool.

		:returns: numpy array (Om0,w0,si8,ns--optionally) 

		"""

		if with_ns:
			return np.array([self.Om0,self.w0,self.sigma8,self.ns])
		else:
			return np.array([self.Om0,self.w0,self.sigma8])


	def squeezePower(self,power=0.5,with_ns=False):

		"""
		Same as squeeze, but returns the doublet (w0,si8 x Om0^power) where power can be specified

		:param power: power at which to raise Om0 in the combination with sigma8
		:type power: float.

		:param with_ns: if True returns also ns as the last parameter
		:type with_ns: bool.

		:returns: numpy array (w0,si8 x Om^power,ns--optionally)

		"""

		if with_ns:
			return np.array([self.w0,self.sigma8*(self.Om0**power),self.ns])
		else:
			return np.array([self.w0,self.sigma8*(self.Om0**power)])

	@classmethod
	def getModels(cls,root_path="/default"):

		"""
		On call, this class method returns a list of IGS1 instances initialized with the cosmological parameters of all the models available in the suite

		:param root_path: path of your IGS1 simulations copy
		:type root_path: str.

		:returns: list of IGS1 instances initialized with the cosmological parameters available in the suite

		"""

		#These are the available models in the suite
		available_models = [ {"Om0":0.26,"Ode0":0.74,"w0":-1.0,"si8":0.798,"ns":0.960} ] 
		available_models.append({"Om0":0.29,"Ode0":0.71,"w0":-1.0,"si8":0.798,"ns":0.960})
		available_models.append({"Om0":0.23,"Ode0":0.77,"w0":-1.0,"si8":0.798,"ns":0.960})
		available_models.append({"Om0":0.26,"Ode0":0.74,"w0":-0.8,"si8":0.798,"ns":0.960})
		available_models.append({"Om0":0.26,"Ode0":0.74,"w0":-1.2,"si8":0.798,"ns":0.960})
		available_models.append({"Om0":0.26,"Ode0":0.74,"w0":-1.0,"si8":0.850,"ns":0.960})
		available_models.append({"Om0":0.26,"Ode0":0.74,"w0":-1.0,"si8":0.750,"ns":0.960})
		available_models.append({"Om0":0.26,"Ode0":0.74,"w0":-1.0,"si8":0.798,"ns":1.000})
		available_models.append({"Om0":0.26,"Ode0":0.74,"w0":-1.0,"si8":0.798,"ns":0.920})

		model_list = [ cls(root_path=root_path,name="Om{0:.3f}_Ol{1:.3f}_w{2:.3f}_ns{3:.3f}_si{4:.3f}".format(model["Om0"],model["Ode0"],model["w0"],model["ns"],model["si8"]),H0=72.0,Om0=model["Om0"],w0=model["w0"],sigma8=model["si8"],ns=model["ns"]) for model in available_models ]

		return model_list



	def getNames(self,realizations,z=1.0,kind="convergence",big_fiducial_set=False):

		"""
		Get the full name of the IGS1 maps, once a redshift, realization identificator and kind are specified

		:param z: redshift plane of the maps, must be one of [1.0,1.5,2.0]
		:type z: float.

		:param realizations: list of realizations to get the names of, the elements must be in [1,1000]
		:type realizations: list. or int.

		:param kind: decide if retrieve convergence or shear maps, must be one of [convergence,shear1,shear2]
		:type kind: str.

		:param big_fiducial_set: set to True if you want to get the names of the bigger fiducial simulation based on 45 N-body simulations
		:type big_fiducial_set: bool.

		"""

		assert type(realizations) == list or type(realizations) == int
		assert z in [1.0,1.5,2.0],"IGS1 doesn't have maps at redshift {0}".format(z)
		assert kind in ["convergence","shear1","shear2"],"You must select one of these: convergence,shear1,shear2"

		if kind=="convergence":
			prefix = "WL-conv"
			direct = "Maps"
		elif kind=="shear1":
			prefix = "Wl-shear1"
			direct = "shear"
		elif kind=="shear2":
			prefix = "Wl-shear2"
			direct = "shear"

		full_path = self._full_path

		if big_fiducial_set:
			assert self.Om0==0.26 and self.w0==-1.0 and self.sigma8==0.798 and self.ns==0.96
			full_path += "_f"

		full_path = os.path.join(full_path,direct)

		if type(realizations) == int:
			return os.path.join(full_path,"{0}_".format(prefix)+self._series_name+"-"+self._box_string+"_"+self._cosmo_id_string+"_"+str(self._lens_plane_size)+"xy_{0:0004d}r_{1}_{2:0004d}z_og.gre.fit".format(realizations,self._plane_id(z),int(z*100)))
		else:
			return [ os.path.join(full_path,"{0}_".format(prefix)+self._series_name+"-"+self._box_string+"_"+self._cosmo_id_string+"_"+str(self._lens_plane_size)+"xy_{0:0004d}r_{1}_{2:0004d}z_og.gre.fit".format(r,self._plane_id(z),int(z*100))) for r in realizations ]


	def load(self,realization,**kwargs):

		"""
		Reads in a specific realization of the convergence field (in FITS format) and returns a ConvergenceMap instance with the loaded map

		:param realization: the specific realization to read
		:type realization: int.

		:param kwargs: the keyword arguments are passed to the getNames method

		:returns: ConvergenceMap instance with the loaded map

		"""

		filename = self.getNames(realization,**kwargs)
		return ConvergenceMap.load(filename,format=self._data_loader) 






