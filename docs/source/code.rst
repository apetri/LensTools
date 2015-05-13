API
***

.. warning::

	This is still in pre-alpha stage, not tested yet on full scale! Use at your own risk!

.. automodule:: lenstools

Convergence maps
================

.. autoclass:: lenstools.image.convergence.ConvergenceMap
	:inherited-members:

.. autoclass:: lenstools.image.convergence.Mask

Shear maps
==========

.. autoclass:: lenstools.image.shear.ShearMap
	:inherited-members:

Statistics
==========

.. autoclass:: lenstools.statistics.ensemble.Ensemble
	:inherited-members: 

Constraining cosmology
======================

.. autoclass:: lenstools.statistics.constraints.Analysis
	:inherited-members:

.. autoclass:: lenstools.statistics.constraints.FisherAnalysis
	:inherited-members:

.. autoclass:: lenstools.statistics.constraints.LikelihoodAnalysis
	:inherited-members:


Confidence contour plotting
===========================

.. autoclass:: lenstools.statistics.contours.ContourPlot
	:inherited-members:

Noise
=====

.. autoclass:: lenstools.image.noise.GaussianNoiseGenerator
	:inherited-members:

Indexing
========

.. automodule:: lenstools.statistics.index

.. autoclass:: lenstools.statistics.index.Indexer

.. autoclass:: lenstools.statistics.index.PowerSpectrum

.. autoclass:: lenstools.statistics.index.Moments

.. autoclass:: lenstools.statistics.index.Peaks

.. autoclass:: lenstools.statistics.index.PDF

.. autoclass:: lenstools.statistics.index.MinkowskiAll

.. autoclass:: lenstools.statistics.index.MinkowskiSingle

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
	:members: path,mkdir,newCollection,getCollection,collections,newTelescopicMapSet,getTelescopicMapSet,telescopicmapsets

.. autoclass:: lenstools.pipeline.simulation.SimulationCollection
	:members: newRealization,getRealization,realizations,newMapSet,getMapSet,writeCAMB,camb2ngenic

.. autoclass:: lenstools.pipeline.simulation.SimulationIC
	:members: newPlaneSet,getPlaneSet,writeNGenIC,writeGadget2

.. autoclass:: lenstools.pipeline.simulation.SimulationPlanes

.. autoclass:: lenstools.pipeline.simulation.SimulationMaps

.. autoclass:: lenstools.pipeline.simulation.SimulationTelescopicMaps

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

.. autoclass:: lenstools.pipeline.settings.TelescopicMapSettings
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


Limber integration
==================

.. autoclass:: lenstools.simulations.limber.LimberIntegrator
	:inherited-members:

Defaults
========

.. automodule:: lenstools.utils.defaults
	:members: load_fits_default_convergence,load_fits_default_shear,default_callback_loader,convergence_measure_all

MPI
===

.. autoclass:: lenstools.utils.mpi.MPIWhirlPool
	:inherited-members:

Fast Fourier Transforms
=======================

.. autoclass:: lenstools.utils.fft.NUMPYFFTPack
	:inherited-members:
