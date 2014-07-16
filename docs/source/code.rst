API
***

.. warning::

	This is still in pre-alpha stage, not tested yet on full scale! Use at your own risk!

.. automodule:: lenstools

Convergence maps
================

.. autoclass:: lenstools.ConvergenceMap
	:inherited-members:

Shear maps
==========

.. autoclass:: lenstools.ShearMap
	:inherited-members:

Statistics
==========

.. autoclass:: lenstools.Ensemble
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

.. autoclass:: lenstools.index.Peaks

.. autoclass:: lenstools.index.MinkowskiAll

.. autoclass:: lenstools.index.MinkowskiSingle

Defaults
========

.. automodule:: lenstools.defaults
	:members: load_fits_default_convergence,load_fits_default_shear,default_callback_loader,convergence_measure_all

Limber integration
==================

.. automodule:: lenstools.limber
	:members: LimberIntegrator