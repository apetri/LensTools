import os
from ..pipeline.simulation import SimulationBatch

from ..pipeline.settings import *
from ..pipeline.remote import LocalGit

from ..pipeline.simulation import LensToolsCosmology
from ..pipeline.settings import Gadget2Settings

from .. import data,dataExtern

import numpy as np

try:
	os.mkdir("SimTest")
except:
	pass

home = "SimTest/Home"
storage = "SimTest/Storage"

cosmologies = [LensToolsCosmology(),LensToolsCosmology(Om0=0.4)]
box_sizes = [15.0,240.0]
part = [32,512]
seeds = [0,11,222]

#Check environment settings
env = EnvironmentSettings(home=home,storage=storage)

#Instantiate a batch of simulations
batch = SimulationBatch(env)

def test_directory_tree():

	#Create two simulation models
	for cosmo in cosmologies:
	
		simulation_model = batch.newModel(cosmology=cosmo,parameters=["Om","Ol","w","si","ns"])

		#Create two different collections, three initial conditions per collection
		for i,box_size in enumerate(box_sizes):
			
			collection = simulation_model.newCollection(box_size=box_size*simulation_model.Mpc_over_h,nside=part[i])

			#Create dummy ngenic power spectrum file
			fp = open(os.path.join(collection.home_subdir,"ngenic_matterpower_z{0:.6f}.txt".format(0.0)),"w")
			fp.close()

			map_settings = MapSettings()
			map_set = collection.newMapSet(map_settings)

			catalog_settings = CatalogSettings()
			catalog = collection.newCatalog(catalog_settings)
			
			for seed in seeds:
				
				ic = collection.newRealization(seed)
				pln_settings = PlaneSettings()
				pln = ic.newPlaneSet(pln_settings)


def test_copy():
	
	try:
		os.mkdir("SimTest/Copy")
	except:
		pass

	bCopy = batch.copyTree("SimTest/Copy")

def test_git():

	try:

		try:
			os.mkdir("SimTest/CopyGit")
		except:
			pass

		git = LocalGit()
		bCopy = batch.copyTree("SimTest/CopyGit",syshandler=git)
		bCopy.commit("Initial directory tree")
	except ImportError:
		pass


def test_present():

	for model in batch.available:
		for collection in model.collections:

			mp = collection.getMapSet("Maps")
			print(mp)

			ct = collection.getCatalog("Catalog")
			print(ct)

			for ic in collection.realizations:
				
				print(ic)
				pln = ic.getPlaneSet("Planes")
				print(pln)


def test_NGenICParam():

	#Create a parameter file for all the initial conditions present in the batch
	settings = NGenICSettings(GlassFile=data("dummy_glass_little_endian.dat"))

	for model in batch.available:
		for collection in model.collections:
			for ic in collection.realizations:
				ic.writeNGenIC(settings)


def test_Gadget2Param():

	#Create a parameter file for all the initial conditions present in the batch
	settings = Gadget2Settings()

	for model in batch.available:
		for collection in model.collections:
			for ic in collection.realizations:
				ic.writeGadget2(settings)
