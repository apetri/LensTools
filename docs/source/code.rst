Documentation for the code
**************************

.. automodule:: lenstools

Limber integration
==================

.. automodule:: lenstools.limber
	:members: LimberIntegrator

Convergence maps
================

.. automodule:: lenstools.topology
	:members: ConvergenceMap

.. autoclass:: lenstools.ConvergenceMap
	:members: __add__

Shear maps
==========

.. automodule:: lenstools.shear
	:members: ShearMap

Statistics
==========

.. automodule:: lenstools.statistics
	:members: Ensemble

.. autoclass:: lenstools.Ensemble
	:members: __add__,__mul__

Noise
=====

.. automodule:: lenstools.noise
	:members: GaussianNoiseGenerator

Indexing
========

.. automodule:: lenstools.index

.. autoclass:: lenstools.index.Indexer

.. autoclass:: lenstools.index.PowerSpectrum

.. autoclass:: lenstools.index.Peaks

Defaults
========

.. automodule:: lenstools.defaults
	:members: load_fits_default_convergence,load_fits_default_shear,default_callback_loader,convergence_measure_all