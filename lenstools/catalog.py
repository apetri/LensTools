"""

.. module:: catalog
	:platform: Unix
	:synopsis: This module handles galaxy catalogs and their properties


.. moduleauthor:: Andrea Petri <apetri@phys.columbia.edu>


"""

import numpy as np
from scipy.ndimage.filters import gaussian_filter
import astropy.table as tbl
import astropy.units as u

try:
	import matplotlib.pyplot as plt
	matplotlib = True
except ImportError:
	matplotlib = False

import extern as ext

##########################################################
################Catalog class#############################
##########################################################

class Catalog(tbl.Table):

	"""
	Class handler of a galaxy catalogue, inherits all the functionality from the astropy.table.Table

	"""

	def __init__(self,*args,**kwargs):

		#Call parent constructor
		super(Catalog,self).__init__(*args,**kwargs)

		#Set spatial information
		self.setSpatialInfo()


	########################################################################################

	def setSpatialInfo(self,field_x="x",field_y="y",unit=u.deg):

		"""
		Sets the spatial information in the catalog

		:param field_x: name of the column that contains the x coordinates of the objects
		:type field_x: str.

		:param field_y: name of the column that contains the y coordinates of the objects
		:type field_y: str.

		:param unit: measure unit of the spatial coordinates
		:type unit: astropy.unit
		"""
		
		self._field_x = field_x
		self._field_y = field_y
		self._position_unit = unit


	########################################################################################

	
	def pixelize(self,map_size,npixel,field_quantity=None,origin=np.zeros(2)*u.deg,smooth=None):

		"""
		Constructs a two dimensional square pixelized map version of one of the scalar properties in the catalog by assigning its objects on a grid

		:param map_size: spatial size of the map
		:type map_size: quantity

		:param npixel: number of pixels on a side
		:type npixel: int.

		:param field_quantity: name of the catalog quantity to map; if None, 1 is assumed
		:type field_quantity: str.

		:param origin: two dimensional coordinates of the origin of the map
		:type origin: array with units

		:param smooth: if not None, the map is smoothed with a gaussian filter of scale smooth
		:type smooth: quantity

		:returns: two dimensional scalar array with the pixelized field (pixels with no objects are treated as NaN)
		:rtype: array 

		"""

		#Safety check
		assert map_size.unit.physical_type==self._position_unit.physical_type
		assert len(origin)==2
		assert origin.unit.physical_type==self._position_unit.physical_type
		assert self._field_x in self.columns,"There is no {0} field in the catalog!".format(self._field_x)
		assert self._field_y in self.columns,"There is no {0} field in the catalog!".format(self._field_y)

		if field_quantity is not None:
			assert field_quantity in self.columns,"There is no {0} field in the catalog!".format(field_quantity)

		#Horizontal and vertical positions
		x = self.columns[self._field_x] - origin[0].to(self._position_unit).value
		y = self.columns[self._field_y] - origin[1].to(self._position_unit).value

		if field_quantity is not None:
			scalar = self.columns[field_quantity].astype(np.float)
		else:
			scalar = np.ones(len(x))

		#Make sure x,y,scalar have all the same length
		assert len(x)==len(y) and len(y)==len(scalar)

		#Perform the pixelization
		scalar_map = np.zeros((npixel,npixel))
		scalar_map.fill(np.nan)
		ext._pixelize.grid2d(x,y,scalar,map_size.to(self._position_unit).value,scalar_map)

		#Maybe smooth
		if smooth is not None:
			
			#Smoothing scale in pixel
			assert smooth.unit.physical_type==self._position_unit.physical_type
			smooth_in_pixel = (smooth * npixel / map_size).decompose().value

			#Replace NaN with zeros
			scalar_map[np.isnan(scalar_map)] = 0.0

			#Smooth
			scalar_map = gaussian_filter(scalar_map,sigma=smooth_in_pixel)

		#Return
		return scalar_map



	########################################################################################

	
	def visualize(self,map_size,npixel,field_quantity=None,origin=np.zeros(2)*u.deg,smooth=None,fig=None,ax=None,colorbar=False):

		"""
		Visualize a two dimensional square pixelized map version of one of the scalar properties in the catalog by assigning its objects on a grid (the pixelization is performed using the pixelize routine)

		:param map_size: spatial size of the map
		:type map_size: quantity

		:param npixel: number of pixels on a side
		:type npixel: int.

		:param field_quantity: name of the catalog quantity to map; if None, 1 is assumed
		:type field_quantity: str.

		:param origin: two dimensional coordinates of the origin of the map
		:type origin: array with units

		:param smooth: if not None, the map is smoothed with a gaussian filter of scale smooth
		:type smooth: quantity

		:returns: two dimensional scalar array with the pixelized field (pixels with no objects are treated as NaN)
		:rtype: array 

		"""

		if not matplotlib:
			raise ImportError("matplotlib is not installed, cannot visualize!")

		#Create figure
		if (fig is None) or (ax is None):
			
			self.fig,self.ax = plt.subplots()

		else:

			self.fig = fig
			self.ax = ax

		#Compute the pixelization
		scalar_map = self.pixelize(map_size,npixel,field_quantity,origin,smooth)

		#Plot 
		ax0 = self.ax.imshow(scalar_map.T,origin="lower",interpolation="nearest",extent=[origin[0].to(self._position_unit).value,(origin[0]+map_size).to(self._position_unit).value,origin[1].to(self._position_unit).value,(origin[1]+map_size).to(self._position_unit).value])
		self.ax.set_xlabel(r"$x$({0})".format(self._position_unit.to_string()))
		self.ax.set_ylabel(r"$y$({0})".format(self._position_unit.to_string()))

		#Colorbar
		if colorbar:
			cbar = plt.colorbar(ax0,ax=self.ax)
			if field_quantity is not None:
				cbar.set_label(field_quantity)
			else:
				cbar.set_label(r"$n$")


	########################################################################################


	def savefig(self,filename):

		"""
		Saves the catalog visualization to an external file

		:param filename: name of the file on which to save the map
		:type filename: str.

		"""

		self.fig.savefig(filename)


	########################################################################################

	def write(self,filename,**kwargs):

		self.meta["NGAL"] = len(self)
		self.meta["AUNIT"] = self._position_unit.to_string()
		super(Catalog,self).write(filename,**kwargs)


##########################################################
################ShearCatalog class########################
##########################################################

class ShearCatalog(Catalog):

	"""
	Class handler of a galaxy shear catalog, inherits all the functionality from the Catalog class

	"""

	########################################################################################

	def write(self,filename,**kwargs):

		self.meta["NGAL"] = len(self)
		super(Catalog,self).write(filename,**kwargs)



