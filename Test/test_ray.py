import sys
sys.path.append("..")

from lenstools.simulations.raytracing import RayTracer,PotentialPlane,DeflectionPlane
import numpy as np
from astropy.units import deg,rad

import logging

logging.basicConfig(level=logging.DEBUG)

def test_ray_simple():

	#Instantiate the RayTracer
	tracer = RayTracer()

	#Add the lenses to the system
	for i in range(2,47):
		tracer.addLens(PotentialPlane.load("../lenstools/simulations/planes/snap{0}_potentialPlane0_normal0.fits".format(i)))

	#Rearrange the lenses according to redshift and roll them randomly along the axes
	tracer.reorderLenses()
	tracer.randomRoll()

	#Start a bucket of light rays from these positions
	b = np.linspace(0.0,2.9,512)
	xx,yy = np.meshgrid(b,b)
	pos = np.array([xx,yy]) * deg

	#Trace the rays
	fin = tracer.shoot(pos,z=2.0)

	#Build the deflection plane
	dfl = DeflectionPlane(fin.value-pos.value,angle=tracer.lens[0].side_angle,redshift=tracer.redshift[-1],cosmology=tracer.lens[0].cosmology,unit=pos.unit)

	#Compute shear and convergence
	conv = dfl.convergence()
	shear = dfl.shear()

	#Finally visualize the result
	conv.visualize(colorbar=True)
	conv.savefig("raytraced_convergence.png")
	shear.visualize(colorbar=True)
	shear.savefig("raytraced_shear.png")