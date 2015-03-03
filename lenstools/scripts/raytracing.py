########################################################
############Ray Tracing scripts#########################
########################################################
from __future__ import division,with_statement

import os
import logging
import time
import cPickle
import resource

from lenstools.utils import MPIWhirlPool

from lenstools.convergence import Spin0
from lenstools import ConvergenceMap,ShearMap

from lenstools.simulations.raytracing import RayTracer,PotentialPlane,DeflectionPlane
from lenstools.pipeline.simulation import SimulationModel
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

	#Instantiate the SimulationModel
	model = SimulationModel.getModel(environment,cosmo_id)

	#Scale the box size to the correct units
	nside,box_size = geometry_id.split("b")
	box_size = float(box_size)*model.Mpc_over_h

	#Get the corresponding simulation collection and map batch handlers
	collection = model.getCollection(box_size,nside)
	map_batch = collection.getMapSet(settings.directory_name)

	#Override the settings with the previously pickled ones, if prompted by user
	if settings.override_with_local:

		local_settings_file = os.path.join(map_batch.home_subdir,"settings.p")

		with open(local_settings_file,"r") as settingsfile:
			settings = cPickle.load(settingsfile)
			assert isinstance(settings,MapSettings)

		if (pool is None) or (pool.is_master()):
			logging.warning("Overriding settings with the previously pickled ones at {0}".format(local_settings_file))

	##################################################################
	##################Settings read###################################
	##################################################################

	#Set random seed to generate the realizations
	np.random.seed(settings.seed)

	#Read map angle,redshift and resolution from the settings
	map_angle = settings.map_angle
	source_redshift = settings.source_redshift
	resolution = settings.map_resolution
	nbody_realizations = settings.mix_nbody_realizations
	cut_points = settings.mix_cut_points
	normals = settings.mix_normals
	map_realizations = settings.lens_map_realizations

	#Decide which map realizations this MPI task will take care of (if pool is None, all of them)
	if pool is None:
		first_map_realization = 0
		last_map_realization = map_realizations
		logging.debug("Generating lensing map realizations from {0} to {1}".format(first_map_realization+1,last_map_realization))
	else:
		assert map_realizations%(pool.size+1)==0,"Perfect load-balancing enforced, map_realizations must be a multiple of the number of MPI tasks!"
		realizations_per_task = map_realizations//(pool.size+1)
		first_map_realization = realizations_per_task*pool.rank
		last_map_realization = realizations_per_task*(pool.rank+1)
		logging.debug("Task {0} will generate lensing map realizations from {1} to {2}".format(pool.rank,first_map_realization+1,last_map_realization))

	#Planes will be read from this path
	plane_path = os.path.join(collection.storage_subdir,"ic{0}",settings.plane_set)
	logging.info("Reading planes from {0}".format(plane_path.format("-".join([str(n) for n in nbody_realizations]))))

	#Read how many snapshots are available
	with open(os.path.join(plane_path.format(nbody_realizations[0]),"info.txt"),"r") as infofile:
		num_snapshots = len(infofile.readlines())

	#Construct the randomization matrix that will differentiate between realizations; this needs to have shape map_realizations x num_snapshots x 3 (ic+cut_points+normals)
	randomizer = np.zeros((map_realizations,num_snapshots,3),dtype=np.int)
	randomizer[:,:,0] = np.random.randint(low=0,high=len(nbody_realizations),size=(map_realizations,num_snapshots))
	randomizer[:,:,1] = np.random.randint(low=0,high=len(cut_points),size=(map_realizations,num_snapshots))
	randomizer[:,:,2] = np.random.randint(low=0,high=len(normals),size=(map_realizations,num_snapshots))

	logging.debug("Randomization matrix has shape {0}".format(randomizer.shape))

	#Save path for the maps
	save_path = os.path.join(map_batch.storage_subdir)
	logging.info("Lensing maps will be saved to {0}".format(save_path))

	begin = time.time()

	#We need one of these for cycles for each map random realization
	for r in range(first_map_realization,last_map_realization):

		#Instantiate the RayTracer
		tracer = RayTracer()

		start = time.time()
		last_timestamp = start

		#############################################################
		###############Add the lenses to the system##################
		#############################################################

		#Open the info file to read the lens specifications (assume the info file is the same for all nbody realizations)
		infofile = open(os.path.join(plane_path.format(nbody_realizations[0]),"info.txt"),"r")

		#Read the info file line by line, and decide if we should add the particular lens corresponding to that line or not
		for s in range(num_snapshots):

			#Read the line
			line = infofile.readline().strip("\n")

			#Stop if there is nothing more to read
			if line=="":
				break

			#Split the line in snapshot,distance,redshift
			line = line.split(",")

			snapshot_number = int(line[0].split("=")[1])
		
			distance,unit = line[1].split("=")[1].split(" ")
			if unit=="Mpc/h":
				distance = float(distance)*model.Mpc_over_h
			else:
				distance = float(distance)*getattr(u,"unit")

			lens_redshift = float(line[2].split("=")[1])

			#Add the lens to the system only if its redshift is < than the one of the source
			logging.info("Adding lens at redshift {0}".format(lens_redshift))
			plane_name = os.path.join(plane_path.format(nbody_realizations[randomizer[r,s,0]]),"snap{0}_potentialPlane{1}_normal{2}.fits".format(snapshot_number,cut_points[randomizer[r,s,1]],normals[randomizer[r,s,2]]))
			tracer.addLens((plane_name,distance,lens_redshift))

		#Close the infofile
		infofile.close()

		now = time.time()
		logging.info("Plane specification reading completed in {0:.3f}s".format(now-start))
		last_timestamp = now

		#Rearrange the lenses according to redshift and roll them randomly along the axes
		tracer.reorderLenses()

		now = time.time()
		logging.info("Reordering completed in {0:.3f}s".format(now-last_timestamp))
		last_timestamp = now

		#Start a bucket of light rays from a regular grid of initial positions
		b = np.linspace(0.0,map_angle.value,resolution)
		xx,yy = np.meshgrid(b,b)
		pos = np.array([xx,yy]) * map_angle.unit

		#Trace the ray deflections
		jacobian = tracer.shoot(pos,z=source_redshift,kind="jacobians")

		now = time.time()
		logging.info("Jacobian ray tracing for realization {0} completed in {1:.3f}s".format(r+1,now-last_timestamp))
		last_timestamp = now

		#Compute shear,convergence and omega from the jacobians
		if settings.convergence:
		
			convMap = ConvergenceMap(data=1.0-0.5*(jacobian[0]+jacobian[3]),angle=map_angle)
			savename = os.path.join(save_path,"WLconv_{0:04d}r.{1}".format(r+1,settings.format))
			logging.info("Saving convergence map to {0}".format(savename)) 
			convMap.save(savename)

		##############################################################################################################################
	
		if settings.shear:
		
			shearMap = ShearMap(data=np.array([0.5*(jacobian[3]-jacobian[0]),-0.5*(jacobian[1]+jacobian[2])]),angle=map_angle)
			savename = os.path.join(save_path,"WLshear_{0:04d}r.{1}".format(r+1,settings.format))
			logging.info("Saving shear map to {0}".format(savename))
			shearMap.save(savename) 

		##############################################################################################################################
	
		if settings.omega:
		
			omegaMap = Spin0(data=-0.5*(jacobian[2]-jacobian[1]),angle=map_angle)
			savename = os.path.join(save_path,"WLomega_{0:04d}r.{1}".format(r+1,settings.format))
			logging.info("Saving omega map to {0}".format(savename))
			omegaMap.save(savename)

		now = time.time()
		logging.info("Weak lensing calculations for realization {0} completed in {1:.3f}s".format(r+1,now-last_timestamp))
		logging.info("Memory usage: {0:.3f} GB".format(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/(1024**3)))
	
	#Safety sync barrier
	if pool is not None:
		pool.comm.Barrier()

	if (pool is None) or (pool.is_master()):	
		now = time.time()
		logging.info("Total runtime {0:.3f}s".format(now-begin))

