#######################################################################
################Cut lens planes out of a Gadget2 snapshot##############
#######################################################################

from __future__ import division

import sys,os
import json

from lenstools.simulations.logs import logdriver,logstderr,peakMemory,peakMemoryAll

from lenstools.pipeline.simulation import SimulationBatch
from lenstools.pipeline.settings import PlaneSettings,PlaneLightConeSettings

import lenstools.simulations
import lenstools.simulations.nbody 

from lenstools.simulations import DensityPlane,PotentialPlane
from lenstools.image.convergence import ConvergenceMap,OmegaMap

from lenstools.simulations.raytracing import RayTracer

from lenstools.utils import MPIWhirlPool
from lenstools import configuration

import numpy as np
from astropy.cosmology import z_at_value

#FFT engine
fftengine = lenstools.simulations.nbody.fftengine

################################################################
################Constant time snapshots#########################
################################################################

def cnstTime(pool,batch,settings,batch_id,override):

	#Safety check
	assert isinstance(pool,MPIWhirlPool) or (pool is None)
	assert isinstance(batch,SimulationBatch)
	assert isinstance(settings,PlaneSettings)

	#Split the id into the model,collection and realization parts
	cosmo_id,geometry_id,realization_id = batch_id.split("|")

	#Get a handle on the simulation model
	model = batch.getModel(cosmo_id)

	#Scale the box size to the correct units
	nside,box_size = geometry_id.split("b")
	box_size = float(box_size)*model.Mpc_over_h

	#Get the realization number
	ic = int(realization_id.strip("ic"))

	#Instantiate the appropriate SimulationIC object
	realization = model.getCollection(box_size,nside).getRealization(ic)
	snapshot_path = realization.snapshot_subdir
	
	#Snapshot handling
	if settings.snapshot_handler not in lenstools.simulations.snapshot_allowed:
		if (pool is None) or (pool.is_master()):
			logdriver.error("Snapshot handler {0} not implemented".format(settings.snapshot_handler))
		sys.exit(1)
	
	snapshot_handler = getattr(lenstools.simulations,settings.snapshot_handler)
	SnapshotFileBase = realization.SnapshotFileBase

	#Log to user
	if (pool is None) or (pool.is_master()):
		logdriver.info("Fast Fourier Transform operations handler: {0}".format(fftengine.__class__.__name__))
		logdriver.info("N-body snapshot handler: {0}".format(settings.snapshot_handler))
		logdriver.info("Reading snapshots from {0}".format(os.path.join(snapshot_path,SnapshotFileBase+"*")))

	#Construct also the SimulationPlanes instance, to handle the current plane batch
	plane_batch = realization.getPlaneSet(settings.directory_name)
	save_path = plane_batch.storage_subdir

	#Override settings from the command line
	if override is not None:
		
		override_dict = json.loads(override)
		for key in override_dict.keys():

			if not hasattr(settings,key):
				raise AttributeError("Plane settings have no such attribute: '{0}'".format(key))

			if (pool is None) or (pool.is_master()):
				logdriver.info("Overriding original setting {0}={1} with command line value {0}={2}".format(key,getattr(settings,key),override_dict[key]))

			setattr(settings,key,override_dict[key])


	#Override with pickled options in storage subdir if prompted by user
	if settings.override_with_local:
		
		local_settings_file = os.path.join(plane_batch.home,"settings.p")
		settings = PlaneSettings.read(local_settings_file)
		assert isinstance(settings,PlaneSettings)

		if (pool is None) or (pool.is_master()):
			logdriver.warning("Overriding settings with the previously pickled ones at {0}".format(local_settings_file))


	if (pool is None) or (pool.is_master()):
		
		logdriver.info("Planes will be saved to {0}".format(save_path))
		#Open the info file to save the planes information
		infofile = open(plane_batch.path("info.txt"),"w")

	#Read from PlaneSettings
	plane_resolution = settings.plane_resolution
	
	first_snapshot = settings.first_snapshot
	last_snapshot = settings.last_snapshot
	if first_snapshot is None:
		snapshots = settings.snapshots
	else:
		snapshots = range(first_snapshot,last_snapshot+1)
	
	cut_points = settings.cut_points
	normals = settings.normals
	thickness = settings.thickness
	thickness_resolution = settings.thickness_resolution
	smooth = settings.smooth
	kind = settings.kind

	#Place holder for the lensing density on the plane
	density_projected = np.empty((plane_resolution,)*2,dtype=np.float32)

	#Open a RMA window on the density placeholder
	if pool is not None:
		pool.openWindow(density_projected)
		
		if pool.is_master():
			logdriver.debug("Opened density window of type {0}".format(pool._window_type))

	#Pre--compute multipoles for solving Poisson equation
	lx,ly = np.meshgrid(fftengine.fftfreq(plane_resolution),fftengine.rfftfreq(plane_resolution),indexing="ij")
	l_squared = lx**2 + ly**2
				
	#Avoid dividing by 0
	l_squared[0,0] = 1.0

	#Arguments
	kwargs = {

	"plane_resolution" : plane_resolution,
	"thickness_resolution" : thickness_resolution,
	"smooth" : smooth,
	"kind" : kind,
	"density_placeholder" : density_projected,
	"l_squared" : l_squared

	}

	#Log the initial memory load
	peak_memory_task,peak_memory_all = peakMemory(),peakMemoryAll(pool)
	if (pool is None) or (pool.is_master()):
		logstderr.info("Initial memory usage: {0:.3f} (task), {1[0]:.3f} (all {1[1]} tasks)".format(peak_memory_task,peak_memory_all))


	num_planes_total = len(snapshots)*len(cut_points)*len(normals)
	nplane = 1 

	#Cycle over each snapshot
	for n in snapshots:

		#Log
		if (pool is None) or (pool.is_master()):
			logdriver.info("Waiting for input files from snapshot {0}...".format(n))
		
		#Open the snapshot
		snapshot_filename = realization.path(snapshot_handler.int2root(SnapshotFileBase,n),where="snapshot_subdir")
		if pool is not None:
			logdriver.info("Task {0} reading nbody snapshot from {1}".format(pool.comm.rank,snapshot_filename))

		snap = snapshot_handler.open(snapshot_filename,pool=pool)

		if pool is not None:
			logdriver.debug("Task {0} read nbody snapshot from {1}".format(pool.comm.rank,snapshot_filename))

		#Get the positions of the particles
		if not hasattr(snap,"positions"):
			snap.getPositions(first=snap._first,last=snap._last)

		#Log memory usage
		if (pool is None) or (pool.is_master()):
			logstderr.debug("Read particle positions: peak memory usage {0:.3f} (task)".format(peakMemory()))

		#Close the snapshot file
		snap.close()

		#Update the summary info file
		if (pool is None) or (pool.is_master()):
			infofile.write("s={0},d={1},z={2}\n".format(n,snap.header["comoving_distance"],snap.header["redshift"]))

		#Cut the lens planes
		for cut,pos in enumerate(cut_points):
			for normal in normals:

				if pool is None or pool.is_master():
					logdriver.info("Cutting {0} plane at {1} with normal {2},thickness {3}, of size {4} x {4}".format(kind,pos,normal,thickness,snap.header["box_size"]))

				############################
				#####Do the cutting#########
				############################
				
				plane,resolution,NumPart = snap.cutPlaneGaussianGrid(normal=normal,center=pos,thickness=thickness,left_corner=np.zeros(3)*snap.Mpc_over_h,**kwargs)
				
				#######################################################################################################################################

				#Save the plane
				plane_file = batch.syshandler.map(os.path.join(save_path,settings.name_format.format(n,kind,cut,normal,settings.format)))

				if (pool is None) or (pool.is_master()):
			
					#Wrap the plane in a PotentialPlane object
					if kind=="potential":
						plane_wrap = PotentialPlane(plane.value,angle=snap.header["box_size"],redshift=snap.header["redshift"],comoving_distance=snap.header["comoving_distance"],cosmology=snap.cosmology,num_particles=NumPart,unit=plane.unit)
					elif kind=="density":
						plane_wrap = DensityPlane(plane,angle=snap.header["box_size"],redshift=snap.header["redshift"],comoving_distance=snap.header["comoving_distance"],cosmology=snap.cosmology,num_particles=NumPart)
					else:
						raise NotImplementedError("Plane of kind '{0}' not implemented!".format(kind))

					#Save the result
					logdriver.info("Saving plane to {0}".format(plane_file))
					plane_wrap.save(plane_file)
					logdriver.debug("Saved plane to {0}".format(plane_file))


				#Log peak memory usage
				peak_memory_task,peak_memory_all = peakMemory(),peakMemoryAll(pool)
				if (pool is None) or (pool.is_master()):
					logstderr.info("Plane {0} of {1} completed, peak memory usage: {2:.3f} (task), {3[0]:.3f} (all {3[1]} tasks)".format(nplane,num_planes_total,peak_memory_task,peak_memory_all))

				nplane += 1
			
				#Safety barrier sync
				if pool is not None:
					pool.comm.Barrier()

	#Safety barrier sync
	if pool is not None:
		pool.comm.Barrier()

		#Close the RMA window
		pool.closeWindow()
		
		if pool.is_master():
			logdriver.debug("Closed density window of type {0}".format(pool._window_type))

	#Close the infofile
	if (pool is None) or (pool.is_master()):
		infofile.close()

	if pool is None or pool.is_master():
		logdriver.info("DONE!!")


#######################################################################
################Light cone projected snapshots#########################
#######################################################################

def lightCone(pool,batch,settings,batch_id,override):

	#Safety check
	assert isinstance(pool,MPIWhirlPool) or (pool is None)
	assert isinstance(batch,SimulationBatch)
	assert isinstance(settings,PlaneLightConeSettings)

	#Split the id into the model,collection and realization parts
	cosmo_id,geometry_id,realization_id = batch_id.split("|")

	#Get a handle on the simulation model
	model = batch.getModel(cosmo_id)

	#Scale the box size to the correct units
	nside,box_size = geometry_id.split("b")
	box_size = float(box_size)*model.Mpc_over_h

	#Get the realization number
	ic = int(realization_id.strip("ic"))

	#Instantiate the appropriate SimulationIC object
	realization = model.getCollection(box_size,nside).getRealization(ic)
	snapshot_path = realization.snapshot_subdir
	
	#Snapshot handling
	if settings.snapshot_handler not in lenstools.simulations.snapshot_allowed:
		if (pool is None) or (pool.is_master()):
			logdriver.error("Snapshot handler {0} not implemented".format(settings.snapshot_handler))
		sys.exit(1)
	
	snapshot_handler = getattr(lenstools.simulations,settings.snapshot_handler)
	SnapshotFileBase = realization.SnapshotFileBase

	#Log to user
	if (pool is None) or (pool.is_master()):
		logdriver.info("Fast Fourier Transform operations handler: {0}".format(fftengine.__class__.__name__))
		logdriver.info("N-body snapshot handler: {0}".format(settings.snapshot_handler))
		logdriver.info("Reading snapshots from {0}".format(os.path.join(snapshot_path,SnapshotFileBase+"*")))

	#Construct also the SimulationPlanes instance, to handle the current plane batch
	plane_batch = realization.getPlaneSet(settings.directory_name)
	save_path = plane_batch.storage_subdir

	#Override settings from the command line
	if override is not None:
		
		override_dict = json.loads(override)
		for key in override_dict.keys():

			if not hasattr(settings,key):
				raise AttributeError("Plane settings have no such attribute: '{0}'".format(key))

			if (pool is None) or (pool.is_master()):
				logdriver.info("Overriding original setting {0}={1} with command line value {0}={2}".format(key,getattr(settings,key),override_dict[key]))

			setattr(settings,key,override_dict[key])


	#Override with pickled options in storage subdir if prompted by user
	if settings.override_with_local:
		
		local_settings_file = os.path.join(plane_batch.home,"settings.p")
		settings = PlaneSettings.read(local_settings_file)
		assert isinstance(settings,PlaneSettings)

		if (pool is None) or (pool.is_master()):
			logdriver.warning("Overriding settings with the previously pickled ones at {0}".format(local_settings_file))


	if (pool is None) or (pool.is_master()):
		
		logdriver.info("Planes will be saved to {0}".format(save_path))
		#Open the info file to save the planes information
		infofile = open(plane_batch.path("info.txt"),"w")

	###############################
	####Read from PlaneSettings####
	###############################

	#Maximum redshift and number of discrete lenses
	zmax = settings.zmax
	num_lenses = settings.num_lenses

	#Normal direction
	normals = (settings.normal,)
	
	#Resolution parameters
	plane_resolution = settings.plane_resolution
	thickness_resolution = settings.thickness_resolution
	smooth = settings.smooth

	#Kind (density, potential or born)
	kind = settings.kind

	#Place holder for the lensing density on the plane
	density_projected = np.empty((plane_resolution,)*2,dtype=np.float32)

	################################################################################
	#Compute the discrete comoving distances which will be the center of the lenses#
	################################################################################

	#Maximum comoving distance and lens thickness
	chi_max = plane_batch.cosmology.comoving_distance(zmax)
	thickness = chi_max/num_lenses

	#Positions of the lenses
	chi_start = thickness/2
	chi_end = chi_max - thickness/2
	chi_centers = np.linspace(chi_start.value,chi_end.value,num_lenses)*chi_max.unit

	#Open a RMA window on the density placeholder
	if pool is not None:
		pool.openWindow(density_projected)
		
		if pool.is_master():
			logdriver.debug("Opened density window of type {0}".format(pool._window_type))

	#Pre--compute multipoles for solving Poisson equation
	lx,ly = np.meshgrid(fftengine.fftfreq(plane_resolution),fftengine.rfftfreq(plane_resolution),indexing="ij")
	l_squared = lx**2 + ly**2
				
	#Avoid dividing by 0
	l_squared[0,0] = 1.0

	#Arguments
	kwargs = {

	"plane_resolution" : plane_resolution,
	"thickness_resolution" : thickness_resolution,
	"smooth" : smooth,
	"density_placeholder" : density_projected,
	"l_squared" : l_squared

	}

	#Log the initial memory load
	peak_memory_task,peak_memory_all = peakMemory(),peakMemoryAll(pool)
	if (pool is None) or (pool.is_master()):
		logstderr.info("Initial memory usage: {0:.3f} (task), {1[0]:.3f} (all {1[1]} tasks)".format(peak_memory_task,peak_memory_all))


	num_planes_total = num_lenses
	nplane = 1 

	###################
	#Open the snapshot#
	###################

	snapshot_filename = realization.path(snapshot_handler.int2root(SnapshotFileBase,0),where="snapshot_subdir")
	if pool is not None:
		logdriver.info("Task {0} reading nbody snapshot from {1}".format(pool.comm.rank,snapshot_filename))

	snap = snapshot_handler.open(snapshot_filename,pool=pool)

	if pool is not None:
		logdriver.debug("Task {0} read nbody snapshot from {1}, particles {2}-{3}".format(pool.comm.rank,snapshot_filename,snap._first,snap._last-1))

	#Get the positions of the particles
	if not hasattr(snap,"positions"):
		snap.getPositions(first=snap._first,last=snap._last)

	#Log memory usage
	if (pool is None) or (pool.is_master()):
		logstderr.debug("Read particle positions: peak memory usage {0:.3f} (task)".format(peakMemory()))

	#Close the snapshot file
	snap.close()

	#############################
	#Check weak lensing settings#
	#############################

	if settings.do_lensing:

		if settings.integration_type=="born":

			if kind!="density":
				if (pool is None or pool.is_master()):
					logdriver.info("Using Born approximation, override {0} lens type with density lens type".format(kind))
				kind = "density"

			tracer = RayTracer(lens_type=DensityPlane)

		else:

			if kind!="potential":
				if (pool is None or pool.is_master()):
					logdriver.info("Integration type {0}, potential required, overriding {1} lens type with potential lens type".format(settings.integration_type,kind))
				kind = "potential"

			tracer = RayTracer(lens_type=PotentialPlane)

	###########################################################################################

	#############################
	#Cycle over the lens centers#
	#############################

	for n,center in enumerate(chi_centers):

		#Compute the redshift
		zlens = z_at_value(plane_batch.cosmology.comoving_distance,center)

		#Update the summary info file
		if (pool is None) or (pool.is_master()):
			infofile.write("s={0},d={1},z={2}\n".format(n,center,zlens))

		#Cut the lens planes
		for cut,pos in enumerate((center,)):
			for normal in normals:

				if pool is None or pool.is_master():
					logdriver.info("Cutting {0} projection at {1} with normal {2},thickness {3}, of size {4} x {4} ({5} x {5})".format(kind,pos,normal,thickness,snap.header["box_size"],plane_resolution))

				############################
				#####Do the cutting#########
				############################
				
				plane,resolution,NumPart = snap.cutPlaneGaussianGrid(normal=normal,center=pos,thickness=thickness,left_corner=np.zeros(3)*snap.Mpc_over_h,kind=kind,**kwargs)
				
				#######################################################################################################################################

				#Save the plane
				plane_file = batch.syshandler.map(os.path.join(save_path,settings.name_format.format(n,kind,cut,normal,settings.format)))

				if (pool is None) or (pool.is_master()):
			
					#Wrap the plane in a Plane object
					if kind=="potential":
						plane_wrap = PotentialPlane(plane.value,angle=snap.header["box_size"],redshift=zlens,comoving_distance=center,cosmology=snap.cosmology,num_particles=NumPart,unit=plane.unit)					
					elif kind=="density":
						plane_wrap = DensityPlane(plane,angle=snap.header["box_size"],redshift=zlens,comoving_distance=center,cosmology=snap.cosmology,num_particles=NumPart)
					else:
						raise NotImplementedError("Lens plane of kind '{0}' not implemented!".format(kind))

					#Add the lens to the tracer if we are doing lensing, otherwise save the lens to disk
					if settings.do_lensing:
						logdriver.info("Adding lens of type {0} at redshift {1:.3f}".format(tracer.lens_type.__name__,plane_wrap.redshift))
						tracer.addLens(plane_wrap)
					else:
						logdriver.info("Saving plane to {0}".format(plane_file))
						plane_wrap.save(plane_file)
						logdriver.debug("Saved plane to {0}".format(plane_file))

				#Log peak memory usage
				peak_memory_task,peak_memory_all = peakMemory(),peakMemoryAll(pool)
				if (pool is None) or (pool.is_master()):
					logstderr.info("Plane {0} of {1} completed, peak memory usage: {2:.3f} (task), {3[0]:.3f} (all {3[1]} tasks)".format(nplane,num_planes_total,peak_memory_task,peak_memory_all))

				nplane += 1
			
				#Safety barrier sync
				if pool is not None:
					pool.comm.Barrier()

	#Safety barrier sync
	if pool is not None:
		pool.comm.Barrier()

		#Close the RMA window
		pool.closeWindow()
		
		if pool.is_master():
			logdriver.debug("Closed density window of type {0}".format(pool._window_type))

	#Close the infofile
	if (pool is None) or (pool.is_master()):
		infofile.close()

	##################################
	##########Raytracing/Born#########
	##################################

	#Only master proceeds
	if settings.do_lensing and ((pool is None) or (pool.is_master())):

		#Add a fudge lens beyond the last redshift
		chi_fudge = chi_end + thickness
		z_fudge = z_at_value(snap.cosmology.comoving_distance,chi_fudge)
		tracer.addLens(tracer.lens_type(np.zeros((plane_resolution,)*2),angle=snap.header["box_size"],redshift=z_fudge,comoving_distance=chi_fudge,cosmology=snap.cosmology,num_particles=None))

		#Reorder the lenses
		tracer.reorderLenses()

		#Create the initial grid of ray initial positions
		b = np.linspace(0.,settings.fov.value,settings.fov_resolution)
		pos = np.array(np.meshgrid(b,b))*settings.fov.unit

		if settings.integration_type=="full":
			jacobian = tracer.shoot(pos,z=zmax,kind="jacobians")
			image = ConvergenceMap(data=1.0-0.5*(jacobian[0]+jacobian[3]),angle=settings.fov,redshift=zmax,cosmology=snap.cosmology)

		elif settings.integration_type=="born":
			image = ConvergenceMap(tracer.convergenceBorn(pos,z=zmax),angle=settings.fov,redshift=zmax,cosmology=snap.cosmology)

		elif settings.integration_type=="postBorn2":
			image = ConvergenceMap(tracer.convergencePostBorn2(pos,z=zmax),angle=settings.fov,redshift=zmax,cosmology=snap.cosmology)

		elif settings.integration_type=="postBorn2-ll":
			image = ConvergenceMap(tracer.convergencePostBorn2(pos,z=zmax,include_gp=False),angle=settings.fov,redshift=zmax,cosmology=snap.cosmology)

		elif settings.integration_type=="postBorn2-gp":
			image = ConvergenceMap(tracer.convergencePostBorn2(pos,z=zmax,include_ll=False),angle=settings.fov,redshift=zmax,cosmology=snap.cosmology)

		elif settings.integration_type=="postBorn1+2":
			image = ConvergenceMap(tracer.convergencePostBorn2(pos,z=zmax,include_first_order=True),angle=settings.fov,redshift=zmax,cosmology=snap.cosmology)

		elif settings.integration_type=="omega2":
			image = OmegaMap(tracer.omegaPostBorn2(pos,z=zmax),angle=settings.fov,redshift=zmax,cosmology=snap.cosmology)

		else:
			raise ValueError("Integration type {0} not supported!".format(settings.integration_type))

		#Save the image
		image_file = batch.syshandler.map(os.path.join(save_path,settings.name_format.format(0,settings.integration_type,zmax,settings.normal,settings.format)))
		logdriver.info("Saving {0} image to {1}".format(image.__class__.__name__,image_file))
		image.save(image_file)
		logdriver.debug("Saved image to {0}".format(image_file))

	#################################################################################################

	#Safety barrier sync
	if pool is not None:
		pool.comm.Barrier()

	#Done
	if (pool is None) or (pool.is_master()):
		logdriver.info("DONE!!")
