"""

.. module:: simulations 
	:platform: Unix
	:synopsis: This module handles the book keeping of simulation sets map names, cosmological parameters, etc... It also provides an API to interact with Gadget2 snapshots and a prototype ray tracing engine


.. moduleauthor:: Andrea Petri <apetri@phys.columbia.edu>


"""

from .design import Design
from .igs1 import IGS1
from .cfhtemu1 import CFHTemu1,CFHTcov
from .gadget2 import Gadget2Snapshot,Gadget2Settings
from .raytracing import Plane,DensityPlane,PotentialPlane,RayTracer
from .nicaea import NicaeaSettings,Nicaea