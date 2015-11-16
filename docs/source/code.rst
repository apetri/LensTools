API
***

.. automodule:: lenstools

Convergence maps
================

.. autoclass:: lenstools.image.convergence.ConvergenceMap
	:inherited-members:

.. autoclass:: lenstools.image.convergence.Mask

Shear maps and catalogs
=======================

.. autoclass:: lenstools.image.shear.ShearMap
	:inherited-members:

.. autoclass:: lenstools.catalog.shear.Catalog
	:members: setSpatialInfo,pixelize,visualize

.. autoclass:: lenstools.catalog.shear.ShearCatalog
	:members: toMap

Statistics
==========

.. autoclass:: lenstools.statistics.ensemble.Ensemble
	:members: nobs,read,readall,read_sql_table,read_sql_query,compute,save,concat,meshgrid,combine_from_dict,combine_columns,suppress_indices,group,covariance,bootstrap,principalComponents,project,compare,selfChi2,shuffle,imshow 

.. autoclass:: lenstools.statistics.database.Database
	:inherited-members:

.. autoclass:: lenstools.statistics.database.ScoreDatabase
	:inherited-members:

.. automodule:: lenstools.statistics.database
	:members: chi2database

Cosmological parameter estimation
=================================

.. autoclass:: lenstools.statistics.constraints.Analysis
	:members: from_features,feature_names,parameter_names,parameter_set,feature_set,parameters,features,add_models,reparametrize,refeaturize,combine_features,find

Fisher matrix calculations
--------------------------

.. autoclass:: lenstools.statistics.constraints.FisherAnalysis
	:members: set_fiducial,fiducial,variations,check,where,varied,compute_derivatives,derivatives,chi2,fit,classify,parameter_covariance,fisher_matrix,ellipse,confidence_ellipse

Non linear parameter dependence
-------------------------------

.. autoclass:: lenstools.statistics.constraints.Emulator
	:members: set_likelihood,train,predict,chi2,chi2Contributions,likelihood,score,approximate_linear

Posterior samplers
------------------

.. automodule:: lenstools.statistics.samplers
	:members: emcee_sampler


Confidence contour plotting
---------------------------

.. autoclass:: lenstools.statistics.contours.ContourPlot
	:inherited-members:


Noise
=====

.. autoclass:: lenstools.image.noise.GaussianNoiseGenerator
	:inherited-members:

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
	:members: parameters,from_specs,write,visualize,savefig,set_title,diagonalCost,cost,sample


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

.. autoclass:: lenstools.simulations.gadget2.Gadget2SnapshotDE

.. autoclass:: lenstools.simulations.gadget2.Gadget2SnapshotPipe

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
	:inherited-members:

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

.. automodule:: lenstools.pipeline.cluster
	:inherited-members:

Real observation sets
=====================

.. automodule:: lenstools.observations
	:inherited-members:

.. autoclass:: lenstools.observations.CFHTLens
	:inherited-members:


Limber integration
==================

.. autoclass:: lenstools.simulations.limber.LimberIntegrator
	:inherited-members:

Defaults
========

.. automodule:: lenstools.utils.defaults
	:members: load_fits_default_convergence,load_fits_default_shear,measure_power_spectrum

MPI
===

.. autoclass:: lenstools.utils.mpi.MPIWhirlPool
	:inherited-members:

Fast Fourier Transforms
=======================

.. autoclass:: lenstools.utils.fft.NUMPYFFTPack
	:inherited-members:
