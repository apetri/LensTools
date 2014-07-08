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
	:members: load_fits_default

.. autoclass:: lenstools.ConvergenceMap
	:members: fromfilename,gradient,hessian,peakCount,minkowskiFunctionals,moments,powerSpectrum

Shear maps
==========

.. automodule:: lenstools.shear
	:members: load_fits_default

.. autoclass:: lenstools.ShearMap
	:members: fromfilename,fromEBmodes,decompose,sticks,convergence,visualizeComponents

Statistics
==========

.. automodule:: lenstools.statistics
	:members: default_callback_loader

.. autoclass:: lenstools.Ensemble
	:members: fromfilelist