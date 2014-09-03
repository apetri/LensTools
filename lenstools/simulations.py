"""

.. module:: simulations 
	:platform: Unix
	:synopsis: This module handles the book keeping of simulation sets map names, cosmological parameters, etc...


.. moduleauthor:: Andrea Petri <apetri@phys.columbia.edu>


"""

from __future__ import division,print_function,with_statement

import os,re

import numpy as np
from astropy.cosmology import FlatwCDM

try:
	from external import _design
	_design = _design
except ImportError:
	_design = None

try:
	import matplotlib.pyplot as plt
	from mpl_toolkits.mplot3d import Axes3D
	matplotlib = True
except ImportError:
	matplotlib = False

######################################
#########Design class#################
######################################

class Design(object):

	"""
	
	A class that proves useful in designing simulation sets: the main functionality provided is the uniform sampling of an arbirtarily high dimensional parameter space. The points in parameter space are chosen to be as spread as possible by minimizing a cost function, but enforcing a latin hypercube structure, in which each parameter value appears once

	"""

	def __init__(self):

		#Input sanity check
		assert _design is not None,"The _design external library has not been installed!"

		#Initialize with 0 points in 0 dimensions
		self.ndim = 0
		self.npoints = 0

		#Useful dictionary containers
		self.parameters = list()
		self.min = dict()
		self.max = dict()
		self.label = dict()
		self.axis = dict()

	def __repr__(self):
		
		if not self.ndim:
			return "This is an empty design!"
		else:
			return "This design has {0} points distributed in a {1}-dimensional parameter space".format(self.npoints,self.ndim)

	def add_parameter(self,parameter_name,min,max,label):

		"""
		Add a dimension to the design by specifying a parameter name, a range and a parameter label (can be in tex format)

		:param parameter_name: the name of the parameter
		:type parameter_name: str.

		:param min: the lower range of the sample interval
		:type min: float.

		:param max: the higher range of the sample interval
		:type max: float.

		:param label: the parameter label you want displayed on a plot, can be in tex format
		:type label: str.

		"""

		assert min<max
		assert parameter_name not in self.parameters,"The parameter is already present!"

		#Fill in containers with the new information
		self.parameters.append(parameter_name)
		self.min[parameter_name] = min
		self.max[parameter_name] = max
		self.label[parameter_name] = label

		#Increase parameter count
		self.axis[parameter_name] = self.ndim
		self.ndim += 1

		#Log the operation
		print("Added a parameter: {0} -> min={1} max={2}".format(parameter_name,min,max))

	def scale(self):
		
		"""
		Scales the points in the design to their respective parameter ranges

		"""

		assert hasattr(self,"points_raw")
		if not hasattr(self,"points"):
			self.points = np.zeros((self.npoints,self.ndim))

		for parameter in self.parameters:
			self.points[:,self.axis[parameter]] = self.min[parameter] + self.points_raw[:,self.axis[parameter]]*(self.max[parameter] - self.min[parameter])


	def put_points(self,npoints):

		"""
		Lay down a number of points on the empty Design: the points are initially layed down on the diagonal of the hypercube

		:param npoints: the number of points to lay down
		:type npoints: int.

		"""
		assert self.ndim>1,"The design must have at least 2 dimensions before laying down points!"
		assert npoints>2, "You must lay down at least 3 points!"

		self.npoints = npoints

		#Lay down points along the diagonal
		self.points_raw = np.outer(np.arange(self.npoints),np.ones(self.ndim)) / (self.npoints - 1)
		
		#Scale to parameter ranges
		self.scale()


	def visualize(self,fig=None,ax=None,parameters=None,**kwargs):

		"""
		Visualize the design configuration using matplotlib

		:param parameters: the parameters to visualize, you can specify two or three of them, by their names. If None, all parameters are visualized
		:type parameters: list.

		"""

		if not matplotlib:
			raise ImportError("matplotlib is not installed, please install it!")

		if parameters is None:
			parameters = self.parameters

		#Check that there are points to plot
		assert hasattr(self,"points"),"There are no points to plot!!"

		#Check that the parameters exist
		for p in parameters:
			assert p in self.parameters,"Parameter {0} is not in your design!".format(p)

		#Check that the parameters to visualize are 2 or 3
		assert len(parameters) in [2,3],"Can plot 2D or 3D projections only!"

		#Instantiate figure and ax objects
		if (fig is None) or (ax is None):
			
			if len(parameters)==2:
				self.fig,self.ax = plt.subplots()
			else:
				self.fig = plt.figure()
				self.ax = self.fig.add_subplot(111,projection="3d")
		
		else:
			
			self.fig,self.ax = fig,ax

		#Lay down the points on the figure
		points = tuple([ self.points[:,self.axis[p]] for p in parameters ])
		self.ax.scatter(*points,**kwargs)

		#Set the labels on the axes
		self.ax.set_xlabel(self.label[parameters[0]])
		self.ax.set_ylabel(self.label[parameters[1]])
		self.ax.set_xlim(self.min[parameters[0]],self.max[parameters[0]])
		self.ax.set_ylim(self.min[parameters[1]],self.max[parameters[1]])

		if len(parameters)==3:
			self.ax.set_zlabel(self.label[parameters[2]])
			self.ax.set_zlim(self.min[parameters[2]],self.max[parameters[2]])

	def savefig(self,filename):

		"""
		Save the visualization to an external file

		:param filename: name of the file on which to save the plot
		:type filename: str.

		"""

		self.fig.savefig(filename)


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


	def __init__(self,H0=70.0,Om0=0.26,w0=-1.0,sigma8=0.798,ns=0.960,root_path=None,name=None):

		super(IGS1,self).__init__(H0,Om0,w0=w0,name=name)
		self.sigma8 = sigma8
		self.ns = ns

		assert root_path is not None,"You must specify the root path of your {0} local copy!".format(self._class_name)

		self.root_path = root_path

		#Don't touch these! 
		self._cosmo_id_string =  "Om{0:.3f}_Ol{1:.3f}_w{2:.3f}_ns{3:.3f}_si{4:.3f}".format(self.Om0,1.0-self.Om0,self.w0,self.ns,self.sigma8)
		self._box_string = str(self._num_particles)+"b"+str(self._box_size_mpc)
		self._full_path = self.root_path.rstrip("/") + "/"+self._series_name+"-"+self._box_string+"_"+self._cosmo_id_string

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

		full_path += "/{0}".format(direct)

		if type(realizations) == int:
			return full_path + "/{0}_".format(prefix)+self._series_name+"-"+self._box_string+"_"+self._cosmo_id_string+"_"+str(self._lens_plane_size)+"xy_{0:0004d}r_{1}_{2:0004d}z_og.gre.fit".format(realizations,self._plane_id(z),int(z*100))
		else:
			return [full_path + "/{0}_".format(prefix)+self._series_name+"-"+self._box_string+"_"+self._cosmo_id_string+"_"+str(self._lens_plane_size)+"xy_{0:0004d}r_{1}_{2:0004d}z_og.gre.fit".format(r,self._plane_id(z),int(z*100)) for r in realizations]


################################
#######EMU1 class###############
################################

class CFHTemu1(IGS1):

	"""
	Class handler of the weak lensing CFHTemu1 simulations set, inherits from IGS1; this simulation suite contains 91 different cosmological models based on 1 N-body simulation each. Each model has 1000 realizations for each of the 13 CFHT subfields 

	"""

	#Don't touch these! 
	_class_name = "CFHTemu1"
	_series_name = "emu1"
	_num_particles = 512
	_box_size_mpc = 240

	@classmethod
	def getModels(cls,root_path="/default"):
		"""
		This class method uses a dictionary file to read in all the cosmological model parameters and instantiate the corresponding CFHTemu1 objects for each one of them

		:param root_path: path of your CFHT emu1 simulations copy
		:type root_path: str.

		:returns: list of CFHTemu1 instances

		"""

		#Build the complete filename
		dict_filename = os.path.join(os.path.dirname(os.path.abspath(__file__)),"book")
		dict_filename = os.path.join(dict_filename,"CFHTemu1.txt")

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









