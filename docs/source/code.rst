API
***

.. warning::

	This is still in pre-alpha stage, not tested yet on full scale! Use at your own risk!

.. automodule:: lenstools

Convergence maps
================

.. autoclass:: lenstools.ConvergenceMap
	:inherited-members:

.. autoclass:: lenstools.Mask

Shear maps
==========

.. autoclass:: lenstools.ShearMap
	:inherited-members:

Statistics
==========

.. autoclass:: lenstools.Ensemble
	:inherited-members: 

Constraining cosmology
======================

.. autoclass:: lenstools.constraints.Analysis
	:inherited-members:

.. autoclass:: lenstools.constraints.FisherAnalysis
	:inherited-members:

.. autoclass:: lenstools.constraints.LikelihoodAnalysis
	:inherited-members:


Confidence contour plotting
===========================

.. autoclass:: lenstools.contours.ContourPlot
	:inherited-members:

Noise
=====

.. autoclass:: lenstools.GaussianNoiseGenerator
	:inherited-members:

Indexing
========

.. automodule:: lenstools.index

.. autoclass:: lenstools.index.Indexer

.. autoclass:: lenstools.index.PowerSpectrum

.. autoclass:: lenstools.index.Moments

.. autoclass:: lenstools.index.Peaks

.. autoclass:: lenstools.index.PDF

.. autoclass:: lenstools.index.MinkowskiAll

.. autoclass:: lenstools.index.MinkowskiSingle

Existing Simulation suites
==========================

.. automodule:: lenstools.simulations

.. autoclass:: lenstools.simulations.IGS1
	:members: getNames, squeeze, load

.. autoclass:: lenstools.simulations.CFHTemu1
	:members: getNames, squeeze, getModels, load

.. autoclass:: lenstools.simulations.CFHTcov
	:members: getNames,squeeze, getModels, load

Simulation design
=================

.. autoclass:: lenstools.simulations.Design
	:inherited-members:


Nicaea bindings
===============

.. autoclass:: lenstools.simulations.NicaeaSettings
	:members: default,available,knobs

.. autoclass:: lenstools.simulations.Nicaea
	:members: fromCosmology,convergencePowerSpectrum,shearTwoPoint

Gadget2 snapshot handling
==========================

.. autoclass:: lenstools.simulations.Gadget2Snapshot
	:inherited-members:

Ray tracing simulations
=======================

.. autoclass:: lenstools.simulations.Plane
	:members: save,load,randomRoll,toReal,toFourier

.. autoclass:: lenstools.simulations.DensityPlane
	:members: potential

.. autoclass:: lenstools.simulations.PotentialPlane
	:members: deflectionAngles,density,shearMatrix

.. autoclass:: lenstools.simulations.raytracing.DeflectionPlane
	:members: jacobian,convergence,shear,omega

.. autoclass:: lenstools.simulations.raytracing.ShearTensorPlane

.. autoclass:: lenstools.simulations.RayTracer
	:inherited-members:

Weak Lensing Simulation Pipeline
================================

Directory tree handling
-----------------------

.. autoclass:: lenstools.pipeline.SimulationBatch
	:members: current,available,info,newModel,getModel,writeCAMBSubmission,writeNGenICSubmission,writeGadget2Submission,writePlaneSubmission,writeRaySubmission

.. autoclass:: lenstools.pipeline.simulation.SimulationModel
	:members: path,mkdir,newCollection,getCollection,collections

.. autoclass:: lenstools.pipeline.simulation.SimulationCollection
	:members: newRealization,getRealization,realizations,newMapSet,getMapSet,writeCAMB,camb2ngenic

.. autoclass:: lenstools.pipeline.simulation.SimulationIC
	:members: newPlaneSet,getPlaneSet,writeNGenIC,writeGadget2

.. autoclass:: lenstools.pipeline.simulation.SimulationPlanes

.. autoclass:: lenstools.pipeline.simulation.SimulationMaps

.. autoclass:: lenstools.pipeline.simulation.SimulationCatalog


Tunable settings
----------------

.. autoclass:: lenstools.pipeline.settings.EnvironmentSettings
	:inherited-members:

.. autoclass:: lenstools.simulations.camb.CAMBSettings
	:inherited-members:

.. autoclass:: lenstools.pipeline.settings.NGenICSettings
	:inherited-members:

.. autoclass:: lenstools.simulations.Gadget2Settings
	:inherited-members:

.. autoclass:: lenstools.pipeline.settings.PlaneSettings
	:inherited-members:

.. autoclass:: lenstools.pipeline.settings.MapSettings
	:inherited-members:

.. autoclass:: lenstools.pipeline.settings.CatalogSettings
	:inherited-members:

.. autoclass:: lenstools.pipeline.settings.JobSettings
	:inherited-members:

Cluster deployment
------------------

.. autoclass:: lenstools.pipeline.deploy.JobHandler
	:inherited-members:

.. autoclass:: lenstools.pipeline.deploy.Directives
	:inherited-members:

.. autoclass:: lenstools.pipeline.deploy.ClusterSpecs
	:inherited-members:

Cluster specific settings
-------------------------

.. autoclass:: lenstools.pipeline.cluster.StampedeHandler

.. autoclass:: lenstools.pipeline.cluster.EdisonHandler


Observation sets
================

.. automodule:: lenstools.observations

.. autoclass:: lenstools.observations.CFHTLens
	:inherited-members:

Defaults
========

.. automodule:: lenstools.defaults
	:members: load_fits_default_convergence,load_fits_default_shear,default_callback_loader,convergence_measure_all

Limber integration
==================

.. autoclass:: lenstools.simulations.limber.LimberIntegrator
	:inherited-members:

MPI
===

.. autoclass:: lenstools.utils.MPIWhirlPool
	:inherited-members:
