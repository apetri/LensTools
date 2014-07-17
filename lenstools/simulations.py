"""

.. module:: simulations 
	:platform: Unix
	:synopsis: This module handles the book keeping of a simulation set map names, cosmological parameters, etc...


.. moduleauthor:: Andrea Petri <apetri@phys.columbia.edu>


"""

from __future__ import division

from astropy.cosmology import FlatwCDM

######################################
###########IGS1 class#################
######################################

class IGS1(FlatwCDM):

	"""
	Class handler of the IGS1 simulations set, inherits the cosmological parameters from the astropy.cosmology.FlatwCDM class

	"""

	def __init__(self,H0=70.0,Om0=0.26,w0=-1.0,sigma8=0.798,ns=0.960,root_path=None):

		super(IGS1,self).__init__(H0,Om0,w0=w0)
		self.sigma8 = sigma8
		self.ns = ns

		assert root_path is not None,"You must specify the root path of your IGS1 local copy!"

		self.root_path = root_path

	def __repr__(self):

		astropy_string = super(IGS1,self).__repr__()
		pieces = astropy_string.split(",")
		si8_piece = u" sigma8={0}".format(self.sigma8)
		ns_piece = u" ns={0}".format(self.ns)

		return ",".join(pieces[:3] + [si8_piece,ns_piece] + pieces[3:])
