from __future__ import division,print_function,with_statement

try:
	from .. import extern as ext
	_design = ext._design
except AttributeError:
	_design = None

from ..statistics.ensemble import Ensemble,Series

import numpy as np
from astropy.table import Table,Column,hstack
import astropy.units as u

try:
	import matplotlib.pyplot as plt
	from mpl_toolkits.mplot3d import Axes3D
	matplotlib = True
except ImportError:
	matplotlib = False

######################################
#########Design class#################
######################################

class Design(Ensemble):

	"""
	
	A class that proves useful in designing simulation sets: the main functionality provided is the uniform sampling of an arbirtarily high dimensional parameter space. The points in parameter space are chosen to be as spread as possible by minimizing a cost function, but enforcing a latin hypercube structure, in which each parameter value appears once

	"""

	##############
	##Properties##
	##############

	@property
	def parameters(self):
		return list(self.columns)

	#########
	###I/O###
	#########

	@classmethod
	def from_specs(cls,npoints,parameters):

		"""
		:param npoints: number of points in the design
		:type npoints: int.

		:param parameters: list of tuples (name,label,min,max)
		:type parameters: list.

		"""

		assert npoints>2

		#Unzip the parameter specs
		par,labels,pmin,pmax = zip(*parameters)
		ndim = len(par)

		#Create the trivial design
		points_raw = np.outer(np.arange(npoints),np.ones(ndim)) / (npoints - 1)
		points = np.zeros_like(points_raw)
		for n,p in enumerate(par):
			points[:,n] = pmin[n] + points_raw[:,n]*(pmax[n] - pmin[n])

		#Instantiate the Design object in the trivial configuration
		trivial_design = cls(points,columns=par)
		
		#Update metadata
		trivial_design._labels = labels
		trivial_design._pmin = pmin
		trivial_design._pmax = pmax
		trivial_design._raw = points_raw
		trivial_design._metadata.extend(["_labels","_pmin","_pmax","_raw"])

		#Return
		return trivial_design


	def write(self,filename=None,max_rows=None,format="ascii.latex",column_format="{0:.3f}",**kwargs):

		"""
		Outputs the points that make up the design in a nicely formatted table

		:param filename: name of the file to which the table will be saved; if None the contents will be printed
		:type filename: str. or file descriptor

		:param max_rows: maximum number of rows in the table, if smaller than the number of points the different chunks are hstacked (useful if there are too many rows for display)
		:type max_rows: int.

		:param format: passed to the Table.write astropy method
		:type format: str.

		:param column_format: format specifier for the numerical values in the Table
		:type column_format: str.

		:param kwargs: the keyword arguments are passed to astropy.Table.write method
		:type kwargs: dict.

		:returns: the Table instance with the design parameters

		"""

		#Check that there is something to save
		assert hasattr(self,"points"),"There are no points in your design yet!"
		names = [ self.label[p] for p in self.parameters ]
		
		if (max_rows is None) or (max_rows>=len(self)):

			#Build the table
			design_table = Table(self.values,names=names)

			#Add the number column to the left
			design_table.add_column(Column(data=range(1,len(self)+1),name=r"$N$"),index=0)

		else:

			#Figure out the splitting
			num_chunks = len(self) // max_rows
			if len(self)%max_rows!=0:
				num_chunks+=1

			#Construct the list of tables to hstack
			design_table = list()

			#Cycle through the chunks and create the sub-tables
			for n in range(num_chunks-1):

				columns = self.values[n*max_rows:(n+1)*max_rows]

				#Build the sub-table
				design_table.append(Table(columns,names=names))

				#Add the number column to the left
				design_table[-1].add_column(Column(data=range(n*max_rows+1,(n+1)*max_rows+1),name=r"$N$"),index=0)

			#Create the last sub-table
			columns = self.values[(num_chunks-1)*max_rows:]
			design_table.append(Table(columns,names=names))
			design_table[-1].add_column(Column(data=range((num_chunks-1)*max_rows+1,len(self)+1),name=r"$N$"),index=0)

			#hstack in a single table
			design_table = hstack(design_table)


		#Tune the format
		for colname in design_table.colnames:
			if not design_table.dtype[colname]==np.int:
				design_table[colname].format = column_format

		#Write the table or return it
		if filename is not None:
			design_table.write(filename,format=format,**kwargs)
			return None
		else:
			return design_table


	def visualize(self,fig=None,ax=None,parameters=None,fontsize=20,**kwargs):

		"""
		Visualize the design configuration using matplotlib

		:param parameters: the parameters to visualize, you can specify two or three of them, by their names. If None, all parameters are visualized
		:type parameters: list.

		"""

		if not matplotlib:
			raise ImportError("matplotlib is not installed, please install it!")

		if parameters is None:
			parameters = self.parameters

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
		points = tuple([ self.values[:,self.parameters.index(p)] for p in parameters ])
		self.ax.scatter(*points,**kwargs)

		#Set the labels on the axes
		px,py = self.parameters.index(parameters[0]),self.parameters.index(parameters[1])
		self.ax.set_xlabel(self._labels[px],fontsize=fontsize)
		self.ax.set_ylabel(self._labels[py],fontsize=fontsize)
		self.ax.set_xlim(self._pmin[px],self._pmax[px])
		self.ax.set_ylim(self._pmin[py],self._pmax[py])

		if len(parameters)==3:
			pz = self.parameters.index(parameters[2])
			self.ax.set_zlabel(self._labels[pz],fontsize=fontsize)
			self.ax.set_zlim(self._pmin[pz],self._pmax[pz])

		#Return
		return self.ax

	def savefig(self,filename):

		"""
		Save the visualization to an external file

		:param filename: name of the file on which to save the plot
		:type filename: str.

		"""

		self.fig.savefig(filename)

	def set_title(self,title):

		"""
		Give a title to the visualized design

		:param title: title of the figure
		:type title: str.

		"""

		self.ax.set_title(title)

	################################################################################################################################
	######################These methods make use of the external _design library to compute and optimize design costs###############
	################################################################################################################################

	def diagonalCost(self,Lambda):

		"""
		Computes the cost function of a diagonal configuration with a specified number of points and a metric parameter lambda; the cost function is calculated on the scaled parameter values, which are always between 0 and 1

		:param Lambda: metric parameter of the cost function; if set to 1.0 the cost function corresponds is the Coulomb potential energy
		:type Lambda: float.

		:returns: the value of the cost function for a diagonal configuration

		"""

		if _design is None:
			raise ImportError("This method requires a working GSL installation!")

		assert len(self)>2,"You must lay down at least 3 points!"
		return _design.diagonalCost(len(self),Lambda)

	def cost(self,p,Lambda):

		"""
		Computes the cost function of the current configuration given the metric parameters (p,Lambda)

		:param Lambda: metric parameter of the cost function; if set to 1.0 the cost function corresponds is the Coulomb potential energy
		:type Lambda: float.

		:param p: metric parameter of the cost function; if set to 2.0 the distances between points are the Euclidean ones
		:type p: float.

		:returns: the value of the cost function

		"""

		if _design is None:
			raise ImportError("This method requires a working GSL installation!")

		assert self.shape[1]>1,"The design must have at least 2 dimensions to lay down points!"
		assert len(self)>2,"You must lay down at least 3 points!"

		return _design.cost(self._raw,p,Lambda)**(1.0/Lambda)

	def sample(self,p=2.0,Lambda=1.0,seed=0,maxIterations=10000):

		"""
		Evenly samples the parameter space by minimizing the cost function computed with the metric parameters (p,Lambda); this operation works inplace

		:param Lambda: metric parameter of the cost function; if set to 1.0 the cost function corresponds is the Coulomb potential energy
		:type Lambda: float.

		:param p: metric parameter of the cost function; if set to 2.0 the distances between points are the Euclidean ones
		:type p: float.

		:param seed: random seed with which the sampler random generator will be initialized
		:type seed: int.

		:param maxIterations: maximum number of iterations that the sampler can perform before stopping
		:type maxIterations: int.

		:returns: the relative change of the cost function the last time it varied during the sampling

		"""

		if _design is None:
			raise ImportError("This method requires a working GSL installation!")

		assert self.shape[1]>1,"The design must have at least 2 dimensions to lay down points!"
		assert len(self)>2,"You must lay down at least 3 points!"

		#Create array that holds the values of the cost function
		self.cost_values = np.ones(maxIterations) * -1.0

		deltaPerc = _design.sample(self._raw,p,Lambda,maxIterations,seed,self.cost_values)
		
		#Scale points to correct units
		points = np.zeros_like(self._raw)
		for n,p in enumerate(self.parameters):
			points[:,n] = self._pmin[n] + self._raw[:,n]*(self._pmax[n] - self._pmin[n])

		self[:] = points

		#Cut the cost_values if we stopped before the maxIterations limit
		cut = (self.cost_values==-1).argmin()
		if cut:
			self.cost_values = self.cost_values[:cut]

		return deltaPerc


	#Sample points in an elliptical parameter space
	@classmethod
	def sample_ellipse(cls,npoints,parameters,center,minor,major,angle,radial_power=0.5,**kwargs):

		"""
		Sample points in a two dimensional parameter space whose boundary has the shape of an ellipse

		:param npoints: number of points in the design
		:type npoints: int.

		:param parameters: name of the parameters, the list must have 2 elements
		:type parameters: list.

		:param center: center of the ellipse
		:type center: tuple.

		:param minor: length of the semi-minor axis
		:type minor: float.

		:param major: length of the semi-major axis
		:type major: float.

		:param angle: rotation angle of the major axis with respect to the horizontal (degrees counterclockwise)
		:type angle: float.

		:param radial_power: this parameter controls the density of points around the center (higher for higher density); 0.5 corresponds to uniform sampling
		:type radial_power: float.

		:param kwargs: the keyword arguments are passed to the sample() method
		:type kwarges: dict.

		:returns:
		:rtype: :py:class:`Design`

		"""

		#Two dimensional rotation matrix
		if hasattr(angle,"unit"):
			angle_rad = angle.to(u.rad)
		else:
			angle_rad = angle*np.pi/180.

		rotator = np.array([[np.cos(angle_rad),np.sin(angle_rad)],[-np.sin(angle_rad),np.cos(angle_rad)]])

		#First sample the parameter space in (r,phi) coordinates
		polar_design = cls.from_specs(npoints,[("r","r",0.,1.),("phi","phi",0.,2*np.pi)])
		polar_design.sample(**kwargs) 

		#Transform into (x,y) coordinates by performing a rotation
		cartesian_design = polar_design.apply(lambda r:Series([major*(r["r"]**radial_power)*np.cos(r["phi"]),minor*(r["r"]**radial_power)*np.sin(r["phi"])],index=parameters),axis=1).dot(Ensemble(rotator,index=parameters,columns=parameters))

		#Shift to the provided center
		for c in range(2):
			cartesian_design[cartesian_design.columns[c]] += center[c]

		#Return to user
		return cartesian_design


