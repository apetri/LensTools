########################################################
############Ray Tracing scripts#########################
########################################################

import os
import logging
import time

from lenstools.utils import MPIWhirlPool

from lenstools.convergence import Spin0
from lenstools import ConvergenceMap,ShearMap

from lenstools.simulations.raytracing import RayTracer,PotentialPlane,DeflectionPlane
from lenstools.pipeline.settings import EnvironmentSettings,MapSettings

import numpy as np
import astropy.units as u


def singleRedshift(pool,environment,settings,id):

	#Safety check
	assert isinstance(pool,MPIWhirlPool) or (pool is None)
	assert isinstance(environment,EnvironmentSettings)
	assert isinstance(settings,MapSettings)

	#Separate the id into cosmo_id and geometry_id
	cosmo_id,geometry_id = id.split("|")

	#Set random seed to generate the realizations
	np.random.seed(settings.seed)

	#Read map angle,redshift and resolution from the settings
	map_angle = settings.map_angle
	redshift = settings.source_redshift
	resolution = settings.map_resolution

	#Instantiate the RayTracer
	tracer = RayTracer()

	start = time.time()
	last_timestamp = start

	#TODO: for now we are mixing only one N-body realization, fix in the future
	plane_path = os.path.join(environment.storage,cosmo_id,geometry_id,"ic{0}".format(settings.mix_nbody_realizations[0]),settings.plane_set)
	logging.info("Reading planes from {0}".format(plane_path))

	#Save path for the maps
	save_path = os.path.join(environment.storage,cosmo_id,geometry_id,settings.directory_name)
	logging.info("Lensing maps will be saved to {0}".format(save_path))


	#Add the lenses to the system 

	#TODO: detect automatically which planes to load
	for i in range(11,59):
	
		plane_name = os.path.join(plane_path,"snap{0}_potentialPlane{1}_normal{2}.fits".format(i,0,0))
		logging.info("Reading plane from {0}...".format(plane_name))
		tracer.addLens(PotentialPlane.load(plane_name))

	now = time.time()
	logging.info("Plane loading completed in {0:.3f}s".format(now-start))
	last_timestamp = now

	#Rearrange the lenses according to redshift and roll them randomly along the axes
	tracer.reorderLenses()

	now = time.time()
	logging.info("Reordering completed in {0:.3f}s".format(now-last_timestamp))
	last_timestamp = now

	tracer.randomRoll()

	now = time.time()
	logging.info("Rolling completed in {0:.3f}s".format(now-last_timestamp))
	last_timestamp = now

	#Start a bucket of light rays from a regular grid of initial positions
	b = np.linspace(0.0,map_angle.value,resolution)
	xx,yy = np.meshgrid(b,b)
	pos = np.array([xx,yy]) * map_angle.unit

	#Trace the ray deflections
	jacobian = tracer.shoot(pos,z=redshift,kind="jacobians")

	now = time.time()
	logging.info("Jacobian ray tracing completed in {0:.3f}s".format(now-last_timestamp))
	last_timestamp = now

	#Compute shear,convergence and omega from the jacobians
	if settings.convergence:
		
		conv = ConvergenceMap(data=1.0-0.5*(jacobian[0]+jacobian[3]),angle=map_angle)
		savename = os.path.join(save_path,"WLconv_0001r.{0}".format(settings.format))
		logging.info("Saving convergence map to {0}".format(savename)) 
		conv.save(savename)

	##############################################################################################################################
	
	if settings.shear:
		
		shear = ShearMap(data=np.array([0.5*(jacobian[3]-jacobian[0]),-0.5*(jacobian[1]+jacobian[2])]),angle=map_angle)
		savename = os.path.join(save_path,"WLshear_0001r.{0}".format(settings.format))
		logging.info("Saving shear map to {0}".format(savename))
		shear.save(savename) 

	##############################################################################################################################
	
	if settings.omega:
		
		omega = Spin0(data=-0.5*(jacobian[2]-jacobian[1]),angle=map_angle)
		savename = os.path.join(save_path,"WLomega_0001r.{0}".format(settings.format))
		logging.info("Saving omega map to {0}".format(savename))
		omega.save(savename)

	now = time.time()
	logging.info("Weak lensing calculations completed in {0:.3f}s".format(now-last_timestamp))
	logging.info("Total runtime {0:.3f}s".format(now-start))

