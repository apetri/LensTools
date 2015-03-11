"""

.. module:: catalog
	:platform: Unix
	:synopsis: This module handles galaxy catalogs and their properties


.. moduleauthor:: Andrea Petri <apetri@phys.columbia.edu>


"""

import astropy.table as tbl

##########################################################
################Catalog class#############################
##########################################################

class Catalog(tbl.Table):

	"""
	Class handler of a galaxy catalogue, inherits all the functionality from the astropy.table.Table

	"""

	def pixelize(self):
		pass