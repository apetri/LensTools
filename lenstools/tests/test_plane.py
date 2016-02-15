import os

from ..simulations import Gadget2Snapshot
from ..simulations.raytracing import PotentialPlane

from .. import dataExtern

import numpy as np
import astropy.units as u


def test_nfw():

	#Create Gadget2Snapshot
	snap = Gadget2Snapshot()

	#Add only one particle
	snap.setPositions(np.ones((1,3),dtype=np.float32)*120.0*u.Mpc)
	snap.weights = np.ones(1,dtype=np.float32)
	snap.virial_radius = np.array([200.0]) * u.Mpc
	snap.concentration = np.array([1.0])

	#Add the header
	snap.setHeaderInfo(box_size=240.0*u.Mpc)

	#Cut the plane
	p,b,n = snap.cutPlaneGaussianGrid(center=120.0*u.Mpc,thickness=240.0*u.Mpc,plane_resolution=256,thickness_resolution=1,left_corner=np.zeros(3)*u.Mpc)

	#Build a PotentialPlane
	pln = PotentialPlane(p/p.max(),snap.header["box_size"],comoving_distance=snap.header["comoving_distance"],unit=None,num_particles=n)
	pln.visualize(colorbar=True)
	pln.savefig("nfw.png")