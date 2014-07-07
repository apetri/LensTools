"""

.. module:: statistics 
	:platform: Unix
	:synopsis: This module implements a set of statistical operations on ensembles of weak lensing maps (shear/convergence)


.. moduleauthor:: Andrea Petri <apetri@phys.columbia.edu>


"""

from __future__ import division

from topology import ConvergenceMap
from shear import ShearMap

import numpy as np
from astropy.io import fits

##########################################
########Ensemble class####################
##########################################

class Ensemble(object):

	"""
	A class that handles statistical operations on weak lensing maps

	>>> from lenstools.statistics import Ensemble

	"""

	def __init__(self):
		print("Under construction!!")