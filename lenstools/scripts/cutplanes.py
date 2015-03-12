#######################################################################
################Cut lens planes out of a Gadget2 snapshot##############
#######################################################################

from __future__ import division

import os
import logging
import cPickle

from lenstools.pipeline.simulation import SimulationBatch
from lenstools.pipeline.settings import PlaneSettings

from lenstools.simulations import Gadget2Snapshot,PotentialPlane
from lenstools.utils import MPIWhirlPool

import numpy as np


#######################################################
################Main execution#########################
#######################################################


def main(pool,batch,settings,id):

	#Safety check
	assert isinstance(pool,MPIWhirlPool) or (pool is None)
	assert isinstance(batch,SimulationBatch)
	assert isinstance(settings,PlaneSettings)

	#Split the id into the model,collection and realization parts
	cosmo_id,geometry_id,realization_id = id.split("|")

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
	
	#Base name of the snapshot files
	SnapshotFileBase = realization.SnapshotFileBase + "_"

	#Log to user
	if (pool is None) or (pool.is_master()):
		logging.info("Reading snapshots from {0}".format(os.path.join(snapshot_path,SnapshotFileBase+"*")))

	#Construct also the SimulationPlanes instance, to handle the current plane batch
	plane_batch = realization.getPlaneSet(settings.directory_name)
	save_path = plane_batch.storage_subdir

	#Override with pickled options in storage subdir if prompted by user
	if settings.override_with_local:
		
		local_settings_file = os.path.join(plane_batch.home_subdir,"settings.p")
		
		with open(local_settings_file,"r") as settingsfile:
			settings = cPickle.load(settingsfile)
			assert isinstance(settings,PlaneSettings)

		if (pool is None) or (pool.is_master()):
			logging.warning("Overriding settings with the previously pickled ones at {0}".format(local_settings_file))


	if (pool is None) or (pool.is_master()):
		
		logging.info("Planes will be saved to {0}".format(save_path))
		#Open the info file to save the planes information
		infofile = open(os.path.join(save_path,"info.txt"),"w")

	#Read from PlaneSettings
	plane_resolution = settings.plane_resolution
	first_snapshot = settings.first_snapshot
	last_snapshot = settings.last_snapshot
	cut_points = settings.cut_points
	normals = settings.normals
	thickness = settings.thickness

	for n in range(first_snapshot,last_snapshot+1):

		#Open the snapshot
		snap = Gadget2Snapshot.open(os.path.join(snapshot_path,SnapshotFileBase+"{0:03d}".format(n)),pool=pool)

		if pool is not None:
			logging.info("Rank {0} reading snapshot from {1}".format(pool.comm.rank,snap.header["files"][0]))

		#Get the positions of the particles
		snap.getPositions()

		#Update the summary info file
		if (pool is None) or (pool.is_master()):
			infofile.write("s={0},d={1},z={2}\n".format(n,snap.header["comoving_distance"],snap.header["redshift"]))

		#Cut the lens planes
		for cut,pos in enumerate(cut_points):
			for normal in normals:

				if pool is not None and pool.is_master():
					logging.info("Cutting plane at {0} with normal {1},thickness {2}, of size {3} x {3}".format(pos,normal,thickness,snap.header["box_size"]))

				#Do the cutting
				plane,resolution,NumPart = snap.cutPlaneGaussianGrid(normal=normal,center=pos,thickness=thickness,left_corner=np.zeros(3)*snap.Mpc_over_h,plane_resolution=plane_resolution,thickness_resolution=1,smooth=1,kind="potential")
				plane_file = os.path.join(save_path,"snap{0}_potentialPlane{1}_normal{2}.{3}".format(n,cut,normal,settings.format))

				if pool is None or pool.is_master():
			
					#Wrap the plane in a PotentialPlane object
					potential_plane = PotentialPlane(plane.value,angle=snap.header["box_size"],redshift=snap.header["redshift"],comoving_distance=snap.header["comoving_distance"],cosmology=snap.cosmology,num_particles=NumPart,unit=plane.unit)

					#Save the result
					logging.info("Saving plane to {0}".format(plane_file))
					potential_plane.save(plane_file)
			
			
				if pool is not None:
			
					#Safety barrier sync
					pool.comm.Barrier()


		#Close the snapshot
		snap.close()

	#Safety barrier sync
	if pool is not None:
		pool.comm.Barrier()

	#Close the infofile
	if (pool is None) or (pool.is_master()):
		infofile.close()

	if pool is None or pool.is_master():
		logging.info("DONE!!")
