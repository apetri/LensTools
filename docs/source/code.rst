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

Simulation suites
=================

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

Gadget2 snapshot handling
==========================

.. autoclass:: lenstools.simulations.Gadget2Snapshot
	:inherited-members:

Ray tracing simulations
=======================

.. autoclass:: lenstools.simulations.PotentialPlane
	:members: save,load,randomRoll,deflectionAngles,toReal,toFourier,density

.. autoclass:: lenstools.simulations.raytracing.DeflectionPlane
	:members: jacobian,convergence,shear

.. autoclass:: lenstools.simulations.RayTracer
	:inherited-members:


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

.. autoclass:: lenstools.LimberIntegrator
	:inherited-members:
