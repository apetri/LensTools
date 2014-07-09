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
	:members: fromfilename,gradient,hessian,peakCount,minkowskiFunctionals,moments,powerSpectrum

Shear maps
==========

.. automodule:: lenstools.shear

.. autoclass:: lenstools.ShearMap
	:members: fromfilename,fromEBmodes,decompose,sticks,convergence,visualizeComponents

Statistics
==========

.. automodule:: lenstools.statistics

.. autoclass:: lenstools.Ensemble
	:members: fromfilelist,load,mean

Defaults
========

.. automodule:: defaults
	:members: load_fits_default_convergence,load_fits_default_shear,default_callback_loader