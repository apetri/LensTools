"""

.. module:: catalog
	:platform: Unix
	:synopsis: This module handles galaxy catalogs and their properties


.. moduleauthor:: Andrea Petri <apetri@phys.columbia.edu>


"""

import numpy as np
import astropy.table as tbl
import astropy.units as u

##########################################################
################Catalog class#############################
##########################################################

class Catalog(tbl.Table):

	"""
	Class handler of a galaxy catalogue, inherits all the functionality from the astropy.table.Table

	"""

	def __init__(self,*args,**kwargs):

		#Call parent constructor
		super(self.__class__,self).__init__(*args,**kwargs)

		#Set spatial information
		self.setSpatialInfo()


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

	
	def pixelize(self,map_size,npixel,field_quantity,origin=np.zeros(2)*u.deg):

		"""
		Constructs a two dimensional square pixelized map version of the catalog by assigning its objects on a grid

		:param map_size: spatial size of the map
		:type map_size: quantity

		:param npixel: number of pixels on a side
		:type npixel: int.

		:param field_quantity: name of the catalog quantity to map
		:type field_quantity: str.

		:param origin: two dimensional coordinates of the origin of the map
		:type origin: array with units

		:returns: two dimensional scalar array with the pixelized field
		:rtype: array 

		"""

		#Safety check
		assert map_size.unit.physical_type==self._position_unit.physical_type
		assert len(origin)==2
		assert origin.unit.physical_type==self._position_unit.physical_type
		assert self._field_x in self.columns,"There is no {0} field in the catalog!".format(self._field_x)
		assert self._field_y in self.columns,"There is no {0} field in the catalog!".format(self._field_y)
		assert field_quantity in self.columns,"There is no {0} field in the catalog!".format(field_quantity)

		#Horizontal and vertical positions
		x = self.columns[self._field_x] - origin[0].to(self._position_unit).value
		y = self.columns[self._field_y] - origin[1].to(self._position_unit).value
		scalar = self.columns[field_quantity].astype(np.float)

		#TODO: Perform the pixelization

		raise NotImplementedError