import sys
sys.path.append("..")

from lenstools.simulations.raytracing import RayTracer,PotentialPlane
import numpy as np
from astropy.units import deg,rad

import time

tracer = RayTracer()

start = time.time()

for i in range(2,47):
	tracer.addLens(PotentialPlane.load("../lenstools/simulations/planes/snap{0}_potentialPlane0_normal0.fits".format(i)))

tracer.reorderLenses()
tracer.randomRoll()

end = time.time() - start

b = np.linspace(0.0,2.9,512)
xx,yy = np.meshgrid(b,b)
pos = np.array([xx,yy]) * deg

fin = tracer.shoot(pos,z=2.0)