#!/usr/bin/env python

import numpy as np
import astropy.units as u
from astropy.cosmology import z_at_value

from lenstools.pipeline.simulation import SimulationBatch,LensToolsCosmology
from lenstools.pipeline.settings import EnvironmentSettings,NGenICSettings,PlaneSettings,MapSettings
from lenstools.simulations.camb import CAMBSettings
from lenstools.simulations.gadget2 import Gadget2Settings

#Geometrical parameters
zmax = 3.1
box_size_Mpc_over_h = 240.0
nside = 512
lens_thickness_Mpc = 120.0

#Code settings: interchangeable
camb = CAMBSettings()
ngenic = NGenICSettings()
gadget2 = Gadget2Settings()

#Code settings: lenstools multi-lens-plane algorithm
planes = PlaneSettings.read("planes.ini")
maps = MapSettings.read("maps.ini")


#Initialize the simulation batch directory tree
def init_batch(home,storage,models):

	#Safety check
	assert isinstance(models,np.ndarray)
	assert models.shape[1]==2

	#Get SimulationBatch handle
	batch = SimulationBatch(home=home,storage=storage)

	#Cycle over parameters and create one model per parameter
	for Om,si8 in models:
	
		#Lay down directory tree
		cosmo = LensToolsCosmology(Om0=Om,Ode0=1-Om,sigma8=si8)
		model = batch.newModel(cosmo,parameters=["Om","Ol","si"])
		collection = model.newCollection(
			box_size=box_size_Mpc_over_h*model.Mpc_over_h,
			nside=nside)
		r = collection.newRealization(seed)

		#Plane and catalog directories
		pln = r.newPlaneSet(planes)
		mp = collection.newMapSet(maps)


	#CAMB settings
	for model in batch.models:
		collection = model.collections[0]
		collection.writeCAMB(z=np.array([0.0]),settings=camb)

	#Compute comoving distance to maximum redshift for each model
	d = list()
	for model in batch.available:
		d.append(model.cosmology.comoving_distance(zmax))

	#Compute lens spacings
	d = np.array([dv.value for dv in d]) * d[0].unit

	#We want to make sure there are lenses up to the maximum of these distances
	lens_distances = np.arange(
		lens_thickness_Mpc,d.max().to(u.Mpc).value + lens_thickness_Mpc,
		lens_thickness_Mpc) * u.Mpc

	#Compute the lens redshifts in each models and write the N-body simulations parameter files
	for model in batch.models:

		#Compute the redshifts of the Gadget snapshots
		z = np.zeros_like(lens_distances.value)
		for n,dlens in enumerate(lens_distances):
			z[n] = z_at_value(model.cosmology.comoving_distance,dlens)

		#Assgn values to gadget settings
		gadget2.OutputScaleFactor = np.sort(1/(1+z))

		#Write parameter files		
		collection = model.collections[0]

		#Convert camb power spectra into ngenic ones
		collection.camb2ngenic(z=0.0)

		#Write NGenIC and Gadget2 parameter files
		r = collection.realizations[0]
		r.writeNGenIC(ngenic)
		r.writeGadget2(gadget2)