import sys
sys.path.append("..")

from lenstools.simulations.raytracing import RayTracer,PotentialPlane,DeflectionPlane
from lenstools import ConvergenceMap,ShearMap
import numpy as np
import matplotlib.pyplot as plt
from astropy.units import deg,rad,arcmin

import logging
import time

logging.basicConfig(level=logging.DEBUG)

#Instantiate the RayTracer
tracer = RayTracer(lens_mesh_size=512)

now = time.time()
last_timestamp = now

#Add the lenses to the system
for i in range(11,57):
	tracer.addLens(PotentialPlane.load("../lenstools/simulations/planes/snap{0}_potentialPlane0_normal0.fits".format(i)))

now = time.time()
logging.info("Plane loading completed in {0:.3f}s".format(now-last_timestamp))
last_timestamp = now

#Rearrange the lenses according to redshift and roll them randomly along the axes
tracer.reorderLenses()
tracer.randomRoll()

now = time.time()
logging.info("Reordering and rolling completed in {0:.3f}s".format(now-last_timestamp))
last_timestamp = now

def test_ray_simple():

	z_final = 2.0

	start = time.time()
	last_timestamp = start

	#Start a bucket of light rays from these positions
	b = np.linspace(0.0,2.9,512)
	xx,yy = np.meshgrid(b,b)
	pos = np.array([xx,yy]) * deg

	#Trace the rays
	fin = tracer.shoot(pos,z=z_final)

	now = time.time()
	logging.info("Ray tracing completed in {0:.3f}s".format(now-last_timestamp))
	last_timestamp = now

	#Build the deflection plane
	dfl = DeflectionPlane(fin.value-pos.value,angle=tracer.lens[0].side_angle,redshift=tracer.redshift[-1],cosmology=tracer.lens[0].cosmology,unit=pos.unit)

	#Compute shear and convergence
	conv = dfl.convergence()
	shear = dfl.shear()
	omega = dfl.omega()

	now = time.time()
	logging.info("Weak lensing calculations completed in {0:.3f}s".format(now-last_timestamp))
	last_timestamp = now

	#Finally visualize the result
	conv.visualize(colorbar=True)
	conv.savefig("raytraced_convergence.png")
	omega.visualize(colorbar=True)
	omega.savefig("raytraced_omega.png")
	shear.visualize(colorbar=True)
	shear.savefig("raytraced_shear.png")

	#We want to plot the power spectrum of the raytraced maps
	fig,ax = plt.subplots()
	l_edges = np.arange(200.0,10000.0,100.0)
	l,Pl = conv.powerSpectrum(l_edges)
	ax.plot(l,l*(l+1)*Pl/(2.0*np.pi),label="From ray positions")

	#And why not, E and B modes too
	figEB,axEB = plt.subplots()
	l,EEl,BBl,EBl = shear.decompose(l_edges)
	axEB.plot(l,l*(l+1)*EEl/(2.0*np.pi),label="EE From ray positions",color="black")
	axEB.plot(l,l*(l+1)*BBl/(2.0*np.pi),label="BB From ray positions",color="green")
	axEB.plot(l,l*(l+1)*np.abs(EBl)/(2.0*np.pi),label="EB From ray positions",color="blue")

	#Now compute the shear and convergence raytracing the actual jacobians (more expensive computationally cause it computes the jacobian at every step)
	finJ = tracer.shoot(pos,z=z_final,kind="jacobians")
	conv = ConvergenceMap(data=1.0-0.5*(finJ[0]+finJ[3]),angle=conv.side_angle)
	shear = ShearMap(data=np.array([0.5*(finJ[3]-finJ[0]),-0.5*(finJ[1]+finJ[2])]),angle=shear.side_angle)

	now = time.time()
	logging.info("Jacobian ray tracing completed in {0:.3f}s".format(now-last_timestamp))
	last_timestamp = now

	#Finally visualize the result
	conv.visualize(colorbar=True)
	conv.savefig("raytraced_convergence_jacobian.png")
	shear.visualize(colorbar=True)
	shear.savefig("raytraced_shear_jacobian.png")

	#We want to plot the power spectrum of the raytraced maps
	l,Pl = conv.powerSpectrum(l_edges)
	ax.plot(l,l*(l+1)*Pl/(2.0*np.pi),label="From Jacobians")
	ax.set_xlabel(r"$l$")
	ax.set_ylabel(r"$l(l+1)P_l/2\pi$")
	ax.set_xscale("log")
	ax.set_yscale("log")
	ax.legend()
	fig.savefig("raytracing_conv_power.png")

	#And why not, E and B modes too
	axEB.plot(l,l*(l+1)*EEl/(2.0*np.pi),label="EE From jacobians",color="black",linestyle="--")
	axEB.plot(l,l*(l+1)*BBl/(2.0*np.pi),label="BB From jacobians",color="green",linestyle="--")
	axEB.plot(l,l*(l+1)*np.abs(EBl)/(2.0*np.pi),label="EB From jacobians",color="blue",linestyle="--")
	axEB.set_xlabel(r"$l$")
	axEB.set_ylabel(r"$l(l+1)P_l/2\pi$")
	axEB.set_xscale("log")
	axEB.set_yscale("log")
	axEB.legend(loc="lower right",prop={"size":10})
	figEB.savefig("raytracing_shear_power.png")

	now = time.time()
	logging.info("Total runtime {0:.3f}s".format(now-start))



def test_distortion():

	#Figures
	fig,ax = plt.subplots(1,2,figsize=(16,8))

	#load unlensed image from png file
	image_unlensed = plt.imread("Data/lensing/lens.png")[:,:,0]
	pos_original = (np.array(np.where(image_unlensed>0)) * 2.0/image_unlensed.shape[0]) * deg
	pos_original = np.roll(pos_original,1,axis=0)
	pos_original[1] *= -1
	pos_original[1] += 2.0*deg
	
	#Plot the original image
	ax[0].scatter(pos_original[0],pos_original[1])
	ax[0].set_xlabel(r"$x$({0})".format(pos_original.unit.to_string()))
	ax[0].set_ylabel(r"$y$({0})".format(pos_original.unit.to_string()))
	ax[0].set_title("Original")

	#Perform forward ray tracing with grid interpolation to compute the image distortion
	pos_apparent = tracer.shootForward(pos_original,z=2.0)

	#Plot the distorted image
	ax[1].scatter(pos_apparent[0],pos_apparent[1])
	ax[1].set_xlabel(r"$x$({0})".format(pos_original.unit.to_string()))
	ax[1].set_ylabel(r"$y$({0})".format(pos_original.unit.to_string()))
	ax[1].set_title("Lensed")

	fig.savefig("lens_distortion.png")

