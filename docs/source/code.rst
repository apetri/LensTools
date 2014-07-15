Documentation for the code
**************************

.. automodule:: lenstools

Limber integration
==================

.. automodule:: lenstools.limber
	:members: LimberIntegrator

.. autoclass:: lenstools.limber.LimberIntegrator
	:members: computeConvergence, writeCAMBSettings

Convergence maps
================

.. automodule:: lenstools.topology
	:members: ConvergenceMap

.. autoclass:: lenstools.ConvergenceMap
	:members: fromfilename,gradient,hessian,peakCount,minkowskiFunctionals,moments,powerSpectrum,__add__

Shear maps
==========

.. automodule:: lenstools.shear
	:members: ShearMap

.. autoclass:: lenstools.ShearMap
	:members: fromfilename,fromEBmodes,decompose,sticks,convergence,visualizeComponents

Statistics
==========

.. automodule:: lenstools.statistics
	:members: Ensemble

.. autoclass:: lenstools.Ensemble
	:members: fromfilelist,load,mean,compare,__add__,__mul__,split

Noise
=====

.. automodule:: lenstools.noise
	:members: GaussianNoiseGenerator

.. autoclass:: lenstools.GaussianNoiseGenerator
	:members: forMap,getShapeNoise,fromConvPower

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