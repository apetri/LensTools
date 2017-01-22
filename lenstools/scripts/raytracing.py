########################################################
############Ray Tracing scripts#########################
########################################################
from __future__ import division,with_statement

import sys,os
import time
import gc

from operator import add
from functools import reduce

from lenstools.simulations.logs import logdriver,logstderr,peakMemory,peakMemoryAll

from lenstools.utils.mpi import MPIWhirlPool

from lenstools import ConvergenceMap,OmegaMap,ShearMap
from lenstools.catalog import Catalog,ShearCatalog

from lenstools.simulations.raytracing import RayTracer,DensityPlane
from lenstools.pipeline.simulation import SimulationBatch
from lenstools.pipeline.settings import MapSettings,TelescopicMapSettings,CatalogSettings

from lenstools.scripts import integration_types

import numpy as np
import astropy.units as u

#############################################################
#########Spilt realizations in subdirectories################
#############################################################

def _subdirectories(num_realizations,realizations_in_subdir):

	assert not(num_realizations%realizations_in_subdir),"The number of realizations in each subdirectory must be the same!"
	s = list()

	if num_realizations==realizations_in_subdir:
		return s

	for c in range(num_realizations//realizations_in_subdir):
		s.append("{0}-{1}".format(c*realizations_in_subdir+1,(c+1)*realizations_in_subdir))

	return s

#####################################################################################
#######Callback to call during raytracing to save the convergence at every step######
#####################################################################################

def convergence_callback(jacobian,tracer,k,realization,angle,map_batch,settings):
	convMap = ConvergenceMap(data=1.0-0.5*(jacobian[0]+jacobian[3]),angle=angle)
	savename = os.path.join(map_batch.storage_subdir,"WLconv_z{0:.2f}_{1:04d}r.{2}".format(tracer.redshift[k],realization+1,settings.format))
	logdriver.debug("Saving convergence map to {0}".format(savename)) 
	convMap.save(savename)

################################################
#######Single redshift ray tracing##############
################################################

def singleRedshift(pool,batch,settings,batch_id):

	#Safety check
	assert isinstance(pool,MPIWhirlPool) or (pool is None)
	assert isinstance(batch,SimulationBatch)

	parts = batch_id.split("|")

	if len(parts)==2:

		assert isinstance(settings,MapSettings)
	
		#Separate the id into cosmo_id and geometry_id
		cosmo_id,geometry_id = parts

		#Get a handle on the model
		model = batch.getModel(cosmo_id)

		#Get the corresponding simulation collection and map batch handlers
		collection = [model.getCollection(geometry_id)]
		map_batch = collection[0].getMapSet(settings.directory_name)
		cut_redshifts = np.array([0.0])

	elif len(parts)==1:

		assert isinstance(settings,TelescopicMapSettings)

		#Get a handle on the model
		model = batch.getModel(parts[0])

		#Get the corresponding simulation collection and map batch handlers
		map_batch = model.getTelescopicMapSet(settings.directory_name)
		collection = map_batch.mapcollections
		cut_redshifts = map_batch.redshifts

	else:
		
		if (pool is None) or (pool.is_master()):
			logdriver.error("Format error in {0}: too many '|'".format(batch_id))
		sys.exit(1)


	#Override the settings with the previously pickled ones, if prompted by user
	if settings.override_with_local:

		local_settings_file = os.path.join(map_batch.home_subdir,"settings.p")
		settings = MapSettings.read(local_settings_file)
		assert isinstance(settings,MapSettings)

		if (pool is None) or (pool.is_master()):
			logdriver.warning("Overriding settings with the previously pickled ones at {0}".format(local_settings_file))

	##################################################################
	##################Settings read###################################
	##################################################################

	#Read map angle,redshift and resolution from the settings
	map_angle = settings.map_angle
	source_redshift = settings.source_redshift
	resolution = settings.map_resolution

	if len(parts)==2:

		#########################
		#Use a single collection#
		#########################

		#Read the plane set we should use
		plane_set = (settings.plane_set,)

		#Randomization
		nbody_realizations = (settings.mix_nbody_realizations,)
		cut_points = (settings.mix_cut_points,)
		normals = (settings.mix_normals,)
		map_realizations = settings.lens_map_realizations

	elif len(parts)==1:

		#######################
		#####Telescopic########
		#######################

		#Check that we have enough info
		for attr_name in ["plane_set","mix_nbody_realizations","mix_cut_points","mix_normals"]:
			if len(getattr(settings,attr_name))!=len(collection):
				if (pool is None) or (pool.is_master()):
					logdriver.error("You need to specify a setting {0} for each collection!".format(attr_name))
				sys.exit(1)

		#Read the plane set we should use
		plane_set = settings.plane_set

		#Randomization
		nbody_realizations = settings.mix_nbody_realizations
		cut_points = settings.mix_cut_points
		normals = settings.mix_normals
		map_realizations = settings.lens_map_realizations



	#Decide which map realizations this MPI task will take care of (if pool is None, all of them)
	try:
		realization_offset = settings.first_realization - 1
	except AttributeError:
		realization_offset = 0

	if pool is None:
		first_map_realization = 0 + realization_offset
		last_map_realization = map_realizations + realization_offset
		realizations_per_task = map_realizations
		logdriver.debug("Generating lensing map realizations from {0} to {1}".format(first_map_realization+1,last_map_realization))
	else:
		assert map_realizations%(pool.size+1)==0,"Perfect load-balancing enforced, map_realizations must be a multiple of the number of MPI tasks!"
		realizations_per_task = map_realizations//(pool.size+1)
		first_map_realization = realizations_per_task*pool.rank + realization_offset
		last_map_realization = realizations_per_task*(pool.rank+1) + realization_offset
		logdriver.debug("Task {0} will generate lensing map realizations from {1} to {2}".format(pool.rank,first_map_realization+1,last_map_realization))

	#Planes will be read from this path
	plane_path = os.path.join("{0}","ic{1}","{2}")

	if (pool is None) or (pool.is_master()):
		for c,coll in enumerate(collection):
			logdriver.info("Reading planes from {0}".format(plane_path.format(coll.storage_subdir,"-".join([str(n) for n in nbody_realizations[c]]),plane_set[c])))

	#Plane info file is the same for all collections
	if (not hasattr(settings,"plane_info_file")) or (settings.plane_info_file is None):
		info_filename = batch.syshandler.map(os.path.join(plane_path.format(collection[0].storage_subdir,nbody_realizations[0][0],plane_set[0]),"info.txt"))
	else:
		info_filename = settings.plane_info_file

	if (pool is None) or (pool.is_master()):
		logdriver.info("Reading lens plane summary information from {0}".format(info_filename))

	#Read how many snapshots are available
	with open(info_filename,"r") as infofile:
		num_snapshots = len(infofile.readlines())

	#Save path for the maps
	save_path = map_batch.storage_subdir

	if (pool is None) or (pool.is_master()):
		logdriver.info("Lensing maps will be saved to {0}".format(save_path))

	begin = time.time()

	#Log initial memory load
	peak_memory_task,peak_memory_all = peakMemory(),peakMemoryAll(pool)
	if (pool is None) or (pool.is_master()):
		logstderr.info("Initial memory usage: {0:.3f} (task), {1[0]:.3f} (all {1[1]} tasks)".format(peak_memory_task,peak_memory_all))

	#We need one of these for cycles for each map random realization
	for rloc,r in enumerate(range(first_map_realization,last_map_realization)):

		#Set random seed to generate the realizations
		np.random.seed(settings.seed + r)

		#Instantiate the RayTracer
		tracer = RayTracer()

		#Force garbage collection
		gc.collect()

		#Start timestep
		start = time.time()
		last_timestamp = start

		#############################################################
		###############Add the lenses to the system##################
		#############################################################

		#Open the info file to read the lens specifications (assume the info file is the same for all nbody realizations)
		infofile = open(info_filename,"r")

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

			#Select the right collection
			for n,z in enumerate(cut_redshifts):
				if lens_redshift>=z:
					c = n

			#Randomization of planes
			nbody = np.random.randint(low=0,high=len(nbody_realizations[c]))
			cut = np.random.randint(low=0,high=len(cut_points[c]))
			normal = np.random.randint(low=0,high=len(normals[c]))

			#Log to user
			logdriver.debug("Realization,snapshot=({0},{1}) --> NbodyIC,cut_point,normal=({2},{3},{4})".format(r,s,nbody_realizations[c][nbody],cut_points[c][cut],normals[c][normal]))

			#Add the lens to the system
			logdriver.info("Adding lens at redshift {0}".format(lens_redshift))
			plane_name = batch.syshandler.map(os.path.join(plane_path.format(collection[c].storage_subdir,nbody_realizations[c][nbody],plane_set[c]),settings.plane_name_format.format(snapshot_number,cut_points[c][cut],normals[c][normal],settings.plane_format)))
			tracer.addLens((plane_name,distance,lens_redshift))

		#Close the infofile
		infofile.close()

		now = time.time()
		logdriver.info("Plane specification reading completed in {0:.3f}s".format(now-start))
		last_timestamp = now

		#Rearrange the lenses according to redshift and roll them randomly along the axes
		tracer.reorderLenses()

		now = time.time()
		logdriver.info("Reordering completed in {0:.3f}s".format(now-last_timestamp))
		last_timestamp = now

		#Start a bucket of light rays from a regular grid of initial positions
		b = np.linspace(0.0,map_angle.value,resolution)
		xx,yy = np.meshgrid(b,b)
		pos = np.array([xx,yy]) * map_angle.unit

		if settings.tomographic_convergence:

			#Trace the ray deflections and save the convergence at every step
			tracer.shoot(pos,z=source_redshift,kind="jacobians",callback=convergence_callback,realization=r,angle=map_angle,map_batch=map_batch,settings=settings)

		else:

			#Trace the ray deflections
			jacobian = tracer.shoot(pos,z=source_redshift,kind="jacobians")

			now = time.time()
			logdriver.info("Jacobian ray tracing for realization {0} completed in {1:.3f}s".format(r+1,now-last_timestamp))
			last_timestamp = now

			#Compute shear,convergence and omega from the jacobians
			if settings.convergence or settings.reduced_shear or settings.reduced_shear_convergence:
		
				convMap = ConvergenceMap(data=1.0-0.5*(jacobian[0]+jacobian[3]),angle=map_angle,cosmology=map_batch.cosmology,redshift=source_redshift)
				
				if settings.convergence:
					savename = batch.syshandler.map(os.path.join(save_path,"WLconv_z{0:.2f}_{1:04d}r.{2}".format(source_redshift,r+1,settings.format)))
					logdriver.info("Saving convergence map to {0}".format(savename)) 
					convMap.save(savename)
					logdriver.debug("Saved convergence map to {0}".format(savename)) 

			##############################################################################################################################
	
			if settings.shear or settings.convergence_ks or settings.reduced_shear or settings.reduced_shear_convergence:
		
				shearMap = ShearMap(data=np.array([0.5*(jacobian[3]-jacobian[0]),-0.5*(jacobian[1]+jacobian[2])]),angle=map_angle,cosmology=map_batch.cosmology,redshift=source_redshift)

				if settings.shear:
					savename = batch.syshandler.map(os.path.join(save_path,"WLshear_z{0:.2f}_{1:04d}r.{2}".format(source_redshift,r+1,settings.format)))
					logdriver.info("Saving shear map to {0}".format(savename))
					shearMap.save(savename)

				if settings.convergence_ks:
					convMap = shearMap.convergence() 
					savename = batch.syshandler.map(os.path.join(save_path,"WLconv-ks_z{0:.2f}_{1:04d}r.{2}".format(source_redshift,r+1,settings.format)))
					logdriver.info("Saving convergence (KS) map to {0}".format(savename))
					convMap.save(savename)

				if settings.reduced_shear or settings.reduced_shear_convergence:
					for ng in (0,1):
						shearMap.data[ng] /= (1. - convMap.data)
					
					if settings.reduced_shear:
						savename = batch.syshandler.map(os.path.join(save_path,"WLredshear_z{0:.2f}_{1:04d}r.{2}".format(source_redshift,r+1,settings.format)))
						logdriver.info("Saving reduced shear map to {0}".format(savename))
						shearMap.save(savename)

					if settings.reduced_shear_convergence:
						convMap = shearMap.convergence()
						savename = batch.syshandler.map(os.path.join(save_path,"WLredconv_z{0:.2f}_{1:04d}r.{2}".format(source_redshift,r+1,settings.format)))
						logdriver.info("Saving reduced shear corrected convergence map to {0}".format(savename))
						convMap.save(savename)

			##############################################################################################################################
	
			if settings.omega:
		
				omegaMap = OmegaMap(data=-0.5*(jacobian[2]-jacobian[1]),angle=map_angle,cosmology=map_batch.cosmology,redshift=source_redshift)
				savename = batch.syshandler.map(os.path.join(save_path,"WLomega_z{0:.2f}_{1:04d}r.{2}".format(source_redshift,r+1,settings.format)))
				logdriver.info("Saving omega map to {0}".format(savename))
				omegaMap.save(savename)

		now = time.time()
		
		#Log peak memory usage to stdout
		peak_memory_task,peak_memory_all = peakMemory(),peakMemoryAll(pool)
		logdriver.info("Weak lensing calculations for realization {0} completed in {1:.3f}s".format(r+1,now-last_timestamp))
		logdriver.info("Peak memory usage: {0:.3f} (task), {1[0]:.3f} (all {1[1]} tasks)".format(peak_memory_task,peak_memory_all))

		#Log progress and peak memory usage to stderr
		if (pool is None) or (pool.is_master()):
			logstderr.info("Progress: {0:.2f}%, peak memory usage: {1:.3f} (task), {2[0]:.3f} (all {2[1]} tasks)".format(100*(rloc+1.)/realizations_per_task,peak_memory_task,peak_memory_all))
	
	#Safety sync barrier
	if pool is not None:
		pool.comm.Barrier()

	if (pool is None) or (pool.is_master()):	
		now = time.time()
		logdriver.info("Total runtime {0:.3f}s".format(now-begin))

############################################################################################################################################################################

#########################################################
#######Save intermediate results of LOS integration######
#########################################################

def save_intermediate(add_on,tracer,k,ctype,map_batch=None,map_angle=None,realization=None):
	savename = os.path.join(map_batch.storage,"{0}-lens{1}-{2:04d}r.fits".format(ctype,k,realization))
	logdriver.info("Saving z={0:.3f} add-on to convergence to {1}".format(tracer.redshift[k],savename))
	side = map_angle.to(u.rad).value*tracer.distance[k]
	ConvergenceMap(add_on,angle=side).save(savename)

########################################################################################################################

def losIntegrate(pool,batch,settings,batch_id):

	#Safety check
	assert isinstance(pool,MPIWhirlPool) or (pool is None)
	assert isinstance(batch,SimulationBatch)

	parts = batch_id.split("|")

	if len(parts)==2:

		assert isinstance(settings,MapSettings)
	
		#Separate the id into cosmo_id and geometry_id
		cosmo_id,geometry_id = parts

		#Get a handle on the model
		model = batch.getModel(cosmo_id)

		#Get the corresponding simulation collection and map batch handlers
		collection = [model.getCollection(geometry_id)]
		map_batch = collection[0].getMapSet(settings.directory_name)
		cut_redshifts = np.array([0.0])

	elif len(parts)==1:

		assert isinstance(settings,TelescopicMapSettings)

		#Get a handle on the model
		model = batch.getModel(parts[0])

		#Get the corresponding simulation collection and map batch handlers
		map_batch = model.getTelescopicMapSet(settings.directory_name)
		collection = map_batch.mapcollections
		cut_redshifts = map_batch.redshifts

	else:
		
		if (pool is None) or (pool.is_master()):
			logdriver.error("Format error in {0}: too many '|'".format(batch_id))
		sys.exit(1)


	#Override the settings with the previously pickled ones, if prompted by user
	if settings.override_with_local:

		local_settings_file = os.path.join(map_batch.home_subdir,"settings.p")
		settings = MapSettings.read(local_settings_file)
		assert isinstance(settings,MapSettings)

		if (pool is None) or (pool.is_master()):
			logdriver.warning("Overriding settings with the previously pickled ones at {0}".format(local_settings_file))

	##################################################################
	##################Settings read###################################
	##################################################################

	#Read map angle,redshift and resolution from the settings
	map_angle = settings.map_angle
	source_redshift = settings.source_redshift
	resolution = settings.map_resolution

	if len(parts)==2:

		#########################
		#Use a single collection#
		#########################

		#Read the plane set we should use
		plane_set = (settings.plane_set,)

		#Randomization
		nbody_realizations = (settings.mix_nbody_realizations,)
		cut_points = (settings.mix_cut_points,)
		normals = (settings.mix_normals,)
		map_realizations = settings.lens_map_realizations

	elif len(parts)==1:

		#######################
		#####Telescopic########
		#######################

		#Check that we have enough info
		for attr_name in ["plane_set","mix_nbody_realizations","mix_cut_points","mix_normals"]:
			if len(getattr(settings,attr_name))!=len(collection):
				if (pool is None) or (pool.is_master()):
					logdriver.error("You need to specify a setting {0} for each collection!".format(attr_name))
				sys.exit(1)

		#Read the plane set we should use
		plane_set = settings.plane_set

		#Randomization
		nbody_realizations = settings.mix_nbody_realizations
		cut_points = settings.mix_cut_points
		normals = settings.mix_normals
		map_realizations = settings.lens_map_realizations

	#Integration type
	if settings.integration_type not in integration_types:
		if (pool is None) or (pool.is_master()):
			logdriver.error("Integration type {0} not supported, please choose one in {1}".format(settings.integration_type,integration_types))
		sys.exit(1)

	if (pool is None) or (pool.is_master()):
		logdriver.info("Line of sight integration type: {0}".format(settings.integration_type))

	#Decide which map realizations this MPI task will take care of (if pool is None, all of them)
	try:
		realization_offset = settings.first_realization - 1
	except AttributeError:
		realization_offset = 0

	if pool is None:
		first_map_realization = 0 + realization_offset
		last_map_realization = map_realizations + realization_offset
		realizations_per_task = map_realizations
		logdriver.debug("Generating lensing map realizations from {0} to {1}".format(first_map_realization+1,last_map_realization))
	else:
		assert map_realizations%(pool.size+1)==0,"Perfect load-balancing enforced, map_realizations must be a multiple of the number of MPI tasks!"
		realizations_per_task = map_realizations//(pool.size+1)
		first_map_realization = realizations_per_task*pool.rank + realization_offset
		last_map_realization = realizations_per_task*(pool.rank+1) + realization_offset
		logdriver.debug("Task {0} will generate lensing map realizations from {1} to {2}".format(pool.rank,first_map_realization+1,last_map_realization))

	#Planes will be read from this path
	plane_path = os.path.join("{0}","ic{1}","{2}")

	if (pool is None) or (pool.is_master()):
		for c,coll in enumerate(collection):
			logdriver.info("Reading planes from {0}".format(plane_path.format(coll.storage_subdir,"-".join([str(n) for n in nbody_realizations[c]]),plane_set[c])))

	#Plane info file is the same for all collections
	if (not hasattr(settings,"plane_info_file")) or (settings.plane_info_file is None):
		info_filename = batch.syshandler.map(os.path.join(plane_path.format(collection[0].storage_subdir,nbody_realizations[0][0],plane_set[0]),"info.txt"))
	else:
		info_filename = settings.plane_info_file

	if (pool is None) or (pool.is_master()):
		logdriver.info("Reading lens plane summary information from {0}".format(info_filename))

	#Read how many snapshots are available
	with open(info_filename,"r") as infofile:
		num_snapshots = len(infofile.readlines())

	#Save path for the maps
	save_path = map_batch.storage_subdir

	if (pool is None) or (pool.is_master()):
		logdriver.info("Lensing maps will be saved to {0}".format(save_path))

	begin = time.time()

	#Log initial memory load
	peak_memory_task,peak_memory_all = peakMemory(),peakMemoryAll(pool)
	if (pool is None) or (pool.is_master()):
		logstderr.info("Initial memory usage: {0:.3f} (task), {1[0]:.3f} (all {1[1]} tasks)".format(peak_memory_task,peak_memory_all))

	#We need one of these for cycles for each map random realization
	for rloc,r in enumerate(range(first_map_realization,last_map_realization)):

		#Set random seed to generate the realizations
		np.random.seed(settings.seed + r)

		#Instantiate the RayTracer
		if settings.lens_type=="PotentialPlane":
			tracer = RayTracer()
		elif settings.lens_type=="DensityPlane":
			tracer = RayTracer(lens_type=DensityPlane)
		else:
			raise ValueError("Lens type {0} not recognized!".format(settings.lens_type))

		#Force garbage collection
		gc.collect()

		#Start timestep
		start = time.time()
		last_timestamp = start

		#############################################################
		###############Add the lenses to the system##################
		#############################################################

		#Open the info file to read the lens specifications (assume the info file is the same for all nbody realizations)
		infofile = open(info_filename,"r")

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

			#Select the right collection
			for n,z in enumerate(cut_redshifts):
				if lens_redshift>=z:
					c = n

			#Randomization of planes
			nbody = np.random.randint(low=0,high=len(nbody_realizations[c]))
			cut = np.random.randint(low=0,high=len(cut_points[c]))
			normal = np.random.randint(low=0,high=len(normals[c]))

			#Log to user
			logdriver.debug("Realization,snapshot=({0},{1}) --> NbodyIC,cut_point,normal=({2},{3},{4})".format(r,s,nbody_realizations[c][nbody],cut_points[c][cut],normals[c][normal]))

			#Add the lens to the system
			logdriver.info("Adding lens at redshift {0}".format(lens_redshift))
			plane_name = batch.syshandler.map(os.path.join(plane_path.format(collection[c].storage_subdir,nbody_realizations[c][nbody],plane_set[c]),settings.plane_name_format.format(snapshot_number,cut_points[c][cut],normals[c][normal],settings.plane_format)))
			tracer.addLens((plane_name,distance,lens_redshift))

		#Close the infofile
		infofile.close()

		now = time.time()
		logdriver.info("Plane specification reading completed in {0:.3f}s".format(now-start))
		last_timestamp = now

		#Rearrange the lenses according to redshift and roll them randomly along the axes
		tracer.reorderLenses()

		now = time.time()
		logdriver.info("Reordering completed in {0:.3f}s".format(now-last_timestamp))
		last_timestamp = now

		#Start a bucket of light rays from a regular grid of initial positions
		b = np.linspace(0.0,map_angle.value,resolution)
		xx,yy = np.meshgrid(b,b)
		pos = np.array([xx,yy]) * map_angle.unit

		#Save intermediate results
		if settings.tomographic_convergence:
			callback = save_intermediate
		else:
			callback = None

		#Perform the line of sight integration (choose integration type)
		if settings.integration_type=="born":
			image = tracer.convergenceBorn(pos,z=source_redshift,save_intermediate=False)
			img_type = ConvergenceMap

		elif settings.integration_type=="born-rt":
			image = tracer.convergenceBorn(pos,z=source_redshift,real_trajectory=True,save_intermediate=False)
			img_type = ConvergenceMap

		elif settings.integration_type=="postBorn2":
			image = tracer.convergencePostBorn2(pos,z=source_redshift,save_intermediate=False,include_first_order=False,transpose_up_to=settings.transpose_up_to,callback=callback,map_batch=map_batch,map_angle=map_angle,realization=r+1)
			img_type = ConvergenceMap

		elif settings.integration_type=="postBorn2-ll":
			image = tracer.convergencePostBorn2(pos,z=source_redshift,save_intermediate=False,include_first_order=False,include_gp=False,transpose_up_to=settings.transpose_up_to,callback=callback,map_batch=map_batch,map_angle=map_angle,realization=r+1)
			img_type = ConvergenceMap

		elif settings.integration_type=="postBorn2-gp":
			image = tracer.convergencePostBorn2(pos,z=source_redshift,save_intermediate=False,include_first_order=False,include_ll=False,transpose_up_to=settings.transpose_up_to,callback=callback,map_batch=map_batch,map_angle=map_angle,realization=r+1)
			img_type = ConvergenceMap

		elif settings.integration_type=="postBorn1+2":
			image = tracer.convergencePostBorn2(pos,z=source_redshift,save_intermediate=False,include_first_order=True,callback=callback,transpose_up_to=settings.transpose_up_to,map_batch=map_batch,map_angle=map_angle,realization=r+1)
			img_type = ConvergenceMap

		elif settings.integration_type=="postBorn1+2-gp":
			image = tracer.convergencePostBorn2(pos,z=source_redshift,save_intermediate=False,include_first_order=True,include_ll=False,transpose_up_to=settings.transpose_up_to,callback=callback,map_batch=map_batch,map_angle=map_angle,realization=r+1)
			img_type = ConvergenceMap

		elif settings.integration_type=="postBorn1+2-ll":
			image = tracer.convergencePostBorn2(pos,z=source_redshift,save_intermediate=False,include_first_order=True,include_gp=False,transpose_up_to=settings.transpose_up_to,callback=callback,map_batch=map_batch,map_angle=map_angle,realization=r+1)
			img_type = ConvergenceMap			

		elif settings.integration_type=="omega2":
			image = tracer.omegaPostBorn2(pos,z=source_redshift,save_intermediate=False)
			img_type = OmegaMap

		else:
			raise NotImplementedError

		now = time.time()
		logdriver.info("Line of sight integration for realization {0} completed in {1:.3f}s".format(r+1,now-last_timestamp))
		last_timestamp = now

		#Save the image
		savename = batch.syshandler.map(os.path.join(save_path,"{0}_z{1:.2f}_{2:04d}r".format(settings.integration_type,source_redshift,r+1)))
		if settings.transpose_up_to>=0:
			savename += "_t{0}".format(settings.transpose_up_to)
		savename += ".{0}".format(settings.format)

		logdriver.info("Saving {0} map to {1}".format(settings.integration_type,savename))
		img_type(data=image,angle=map_angle,cosmology=map_batch.cosmology,redshift=source_redshift).save(savename)
		logdriver.debug("Saving {0} map to {1}".format(settings.integration_type,savename)) 

		now = time.time()
		
		#Log peak memory usage to stdout
		peak_memory_task,peak_memory_all = peakMemory(),peakMemoryAll(pool)
		logdriver.info("Weak lensing calculations for realization {0} completed in {1:.3f}s".format(r+1,now-last_timestamp))
		logdriver.info("Peak memory usage: {0:.3f} (task), {1[0]:.3f} (all {1[1]} tasks)".format(peak_memory_task,peak_memory_all))

		#Log progress and peak memory usage to stderr
		if (pool is None) or (pool.is_master()):
			logstderr.info("Progress: {0:.2f}%, peak memory usage: {1:.3f} (task), {2[0]:.3f} (all {2[1]} tasks)".format(100*(rloc+1.)/realizations_per_task,peak_memory_task,peak_memory_all))
	
	#Safety sync barrier
	if pool is not None:
		pool.comm.Barrier()

	if (pool is None) or (pool.is_master()):	
		now = time.time()
		logdriver.info("Total runtime {0:.3f}s".format(now-begin))

############################################################################################################################################################################


###############################################
#######Galaxy catalog ray tracing##############
###############################################

def simulatedCatalog(pool,batch,settings,batch_id):

	#Safety check
	assert isinstance(pool,MPIWhirlPool) or (pool is None)
	assert isinstance(batch,SimulationBatch)
	assert isinstance(settings,CatalogSettings)

	#Separate the id into cosmo_id and geometry_id
	cosmo_id,geometry_id = batch_id.split("|")

	#Get a handle on the model
	model = batch.getModel(cosmo_id)

	#Scale the box size to the correct units
	nside,box_size = geometry_id.split("b")
	box_size = float(box_size)*model.Mpc_over_h

	#Get the corresponding simulation collection and catalog handler
	collection = model.getCollection(box_size,nside)
	catalog = collection.getCatalog(settings.directory_name)

	#Override the settings with the previously pickled ones, if prompted by user
	if settings.override_with_local:

		local_settings_file = os.path.join(catalog.home_subdir,"settings.p")
		settings = CatalogSettings.read(local_settings_file)
		assert isinstance(settings,CatalogSettings)

		if (pool is None) or (pool.is_master()):
			logdriver.warning("Overriding settings with the previously pickled ones at {0}".format(local_settings_file))

	##################################################################
	##################Settings read###################################
	##################################################################

	#Read the catalog save path from the settings
	catalog_save_path = catalog.storage_subdir
	if (pool is None) or (pool.is_master()):
		logdriver.info("Lensing catalogs will be saved to {0}".format(catalog_save_path))

	#Read the total number of galaxies to raytrace from the settings
	total_num_galaxies = settings.total_num_galaxies

	#Pre-allocate numpy arrays
	initial_positions = np.zeros((2,total_num_galaxies)) * settings.catalog_angle_unit
	galaxy_redshift = np.zeros(total_num_galaxies)

	#Keep track of the number of galaxies for each catalog
	galaxies_in_catalog = list()

	#Fill in initial positions and redshifts
	for galaxy_position_file in settings.input_files:

		try:
			galaxies_before = reduce(add,galaxies_in_catalog)
		except TypeError:
			galaxies_before = 0
	
		#Read the galaxy positions and redshifts from the position catalog
		if (pool is None) or (pool.is_master()):
			logdriver.info("Reading galaxy positions and redshifts from {0}".format(galaxy_position_file))
		
		position_catalog = Catalog.read(galaxy_position_file)

		if (pool is None) or (pool.is_master()):
			logdriver.info("Galaxy catalog {0} contains {1} galaxies".format(galaxy_position_file,len(position_catalog)))

		#This is just to avoid confusion
		assert position_catalog.meta["AUNIT"]==settings.catalog_angle_unit.to_string(),"Catalog angle units, {0}, do not match with the ones privided in the settings, {1}".format(position_catalog.meta["AUNIT"],settings.catalog_angle_unit.to_string())

		#Keep track of the number of galaxies in the catalog
		galaxies_in_catalog.append(len(position_catalog))

		if (pool is None) or (pool.is_master()):
			#Save a copy of the position catalog to the simulated catalogs directory
			position_catalog.write(os.path.join(catalog_save_path,os.path.basename(galaxy_position_file)),overwrite=True)

		#Fill in initial positions and redshifts
		initial_positions[0,galaxies_before:galaxies_before+len(position_catalog)] = position_catalog["x"] * getattr(u,position_catalog.meta["AUNIT"])
		initial_positions[1,galaxies_before:galaxies_before+len(position_catalog)] = position_catalog["y"] * getattr(u,position_catalog.meta["AUNIT"])
		galaxy_redshift[galaxies_before:galaxies_before+len(position_catalog)] = position_catalog["z"]

	#Make sure that the total number of galaxies matches, and units are correct
	assert reduce(add,galaxies_in_catalog)==total_num_galaxies,"The total number of galaxies in the catalogs, {0}, does not match the number provided in the settings, {1}".format(reduce(add,galaxies_in_catalog),total_num_galaxies)

	##########################################################################################################################################################
	####################################Initial positions and redshifts of galaxies loaded####################################################################
	##########################################################################################################################################################

	#Read the randomization information from the settings
	nbody_realizations = settings.mix_nbody_realizations
	cut_points = settings.mix_cut_points
	normals = settings.mix_normals
	catalog_realizations = settings.lens_catalog_realizations

	if hasattr(settings,"realizations_per_subdirectory"):
		realizations_in_subdir = settings.realizations_per_subdirectory
	else:
		realizations_in_subdir = catalog_realizations

	#Create subdirectories as necessary
	catalog_subdirectory = _subdirectories(catalog_realizations,realizations_in_subdir)
	if (pool is None) or (pool.is_master()):
		
		for d in catalog_subdirectory:
			
			dir_to_make = os.path.join(catalog_save_path,d)
			if not(os.path.exists(dir_to_make)):
				logdriver.info("Creating catalog subdirectory {0}".format(dir_to_make))
				os.mkdir(dir_to_make)

	#Safety barrier sync
	if pool is not None:
		pool.comm.Barrier() 

	#Decide which map realizations this MPI task will take care of (if pool is None, all of them)
	try:
		realization_offset = settings.first_realization - 1
	except AttributeError:
		realization_offset = 0

	if pool is None:
		first_realization = 0 + realization_offset
		last_realization = catalog_realizations + realization_offset
		realizations_per_task = catalog_realizations
		logdriver.debug("Generating lensing catalog realizations from {0} to {1}".format(first_realization+1,last_realization))
	else:
		assert catalog_realizations%(pool.size+1)==0,"Perfect load-balancing enforced, catalog_realizations must be a multiple of the number of MPI tasks!"
		realizations_per_task = catalog_realizations//(pool.size+1)
		first_realization = realizations_per_task*pool.rank + realization_offset
		last_realization = realizations_per_task*(pool.rank+1) + realization_offset
		logdriver.debug("Task {0} will generate lensing catalog realizations from {1} to {2}".format(pool.rank,first_realization+1,last_realization))


	#Planes will be read from this path
	plane_path = os.path.join(collection.storage_subdir,"ic{0}",settings.plane_set)

	if (pool is None) or (pool.is_master()):
		logdriver.info("Reading planes from {0}".format(plane_path.format("-".join([str(n) for n in nbody_realizations]))))

	#Read how many snapshots are available
	with open(batch.syshandler.map(os.path.join(plane_path.format(nbody_realizations[0]),"info.txt")),"r") as infofile:
		num_snapshots = len(infofile.readlines())

	begin = time.time()

	#Log initial memory load
	peak_memory_task,peak_memory_all = peakMemory(),peakMemoryAll(pool)
	if (pool is None) or (pool.is_master()):
		logstderr.info("Initial memory usage: {0:.3f} (task), {1[0]:.3f} (all {1[1]} tasks)".format(peak_memory_task,peak_memory_all))

	#We need one of these for cycles for each map random realization
	for rloc,r in enumerate(range(first_realization,last_realization)):

		#Set random seed to generate the realizations
		np.random.seed(settings.seed + r)

		#Instantiate the RayTracer
		tracer = RayTracer()

		#Force garbage collection
		gc.collect()

		#Start timestep
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

			#Randomization of planes
			nbody = np.random.randint(low=0,high=len(nbody_realizations))
			cut = np.random.randint(low=0,high=len(cut_points))
			normal = np.random.randint(low=0,high=len(normals))

			#Log to user
			logdriver.debug("Realization,snapshot=({0},{1}) --> NbodyIC,cut_point,normal=({2},{3},{4})".format(r,s,nbody_realizations[nbody],cut_points[cut],normals[normal]))

			#Add the lens to the system
			logdriver.info("Adding lens at redshift {0}".format(lens_redshift))
			plane_name = batch.syshandler.map(os.path.join(plane_path.format(nbody_realizations[nbody]),settings.plane_name_format.format(snapshot_number,cut_points[cut],normals[normal],settings.plane_format)))
			tracer.addLens((plane_name,distance,lens_redshift))

		#Close the infofile
		infofile.close()

		now = time.time()
		logdriver.info("Plane specification reading completed in {0:.3f}s".format(now-start))
		last_timestamp = now

		#Rearrange the lenses according to redshift and roll them randomly along the axes
		tracer.reorderLenses()

		now = time.time()
		logdriver.info("Reordering completed in {0:.3f}s".format(now-last_timestamp))
		last_timestamp = now

		#Trace the ray deflections through the lenses
		jacobian = tracer.shoot(initial_positions,z=galaxy_redshift,kind="jacobians")

		now = time.time()
		logdriver.info("Jacobian ray tracing for realization {0} completed in {1:.3f}s".format(r+1,now-last_timestamp))
		last_timestamp = now

		#Build the shear catalog and save it to disk
		if settings.reduced_shear:
			shear_catalog = ShearCatalog([(jacobian[3]-jacobian[0])/(jacobian[3]+jacobian[0]),-(jacobian[1]+jacobian[2])/(jacobian[3]+jacobian[0])],names=("shear1","shear2"))
		else:
			shear_catalog = ShearCatalog([0.5*(jacobian[3]-jacobian[0]),-0.5*(jacobian[1]+jacobian[2])],names=("shear1","shear2"))

		for n,galaxy_position_file in enumerate(settings.input_files):

			try:
				galaxies_before = reduce(add,galaxies_in_catalog[:n])
			except TypeError:
				galaxies_before = 0
		
			#Build savename
			if settings.reduced_shear:
				shear_root = "WLredshear_"
			else:
				shear_root = "WLshear_"

			if len(catalog_subdirectory):
				shear_catalog_savename = batch.syshandler.map(os.path.join(catalog_save_path,catalog_subdirectory[r//realizations_in_subdir],shear_root+os.path.basename(galaxy_position_file.split(".")[0])+"_{0:04d}r.{1}".format(r+1,settings.format)))
			else:
				shear_catalog_savename = batch.syshandler.map(os.path.join(catalog_save_path,shear_root+os.path.basename(galaxy_position_file.split(".")[0])+"_{0:04d}r.{1}".format(r+1,settings.format)))
			
			if settings.reduced_shear:
				logdriver.info("Saving simulated reduced shear catalog to {0}".format(shear_catalog_savename))
			else:
				logdriver.info("Saving simulated shear catalog to {0}".format(shear_catalog_savename))
			
			shear_catalog[galaxies_before:galaxies_before+galaxies_in_catalog[n]].write(shear_catalog_savename,overwrite=True)

		now = time.time()

		#Log peak memory usage to stdout
		peak_memory_task,peak_memory_all = peakMemory(),peakMemoryAll(pool)
		logdriver.info("Weak lensing calculations for realization {0} completed in {1:.3f}s".format(r+1,now-last_timestamp))
		logdriver.info("Peak memory usage: {0:.3f} (task), {1[0]:.3f} (all {1[1]} tasks)".format(peak_memory_task,peak_memory_all))

		#Log progress and peak memory usage to stderr
		if (pool is None) or (pool.is_master()):
			logstderr.info("Progress: {0:.2f}%, peak memory usage: {1:.3f} (task), {2[0]:.3f} (all {2[1]} tasks)".format(100*(rloc+1.)/realizations_per_task,peak_memory_task,peak_memory_all))


	#Safety sync barrier
	if pool is not None:
		pool.comm.Barrier()

	if (pool is None) or (pool.is_master()):	
		now = time.time()
		logdriver.info("Total runtime {0:.3f}s".format(now-begin))






