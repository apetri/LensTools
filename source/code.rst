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

.. autoclass:: lenstools.ConvergenceMap
	:members: fromfilename,gradient,hessian,peakCount,minkowskiFunctionals,moments,powerSpectrum,__add__

Shear maps
==========

.. automodule:: lenstools.shear

.. autoclass:: lenstools.ShearMap
	:members: fromfilename,fromEBmodes,decompose,sticks,convergence,visualizeComponents

Statistics
==========

.. automodule:: lenstools.statistics

.. autoclass:: lenstools.Ensemble
	:members: fromfilelist,load,mean,__sub__

Noise
=====

.. automodule:: lenstools.noise

.. autoclass:: lenstools.GaussianNoiseGenerator
	:members: forMap,getShapeNoise,fromConvPower

Defaults
========

.. automodule:: lenstools.defaults
	:members: load_fits_default_convergence,load_fits_default_shear,default_callback_loader