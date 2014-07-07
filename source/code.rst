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

.. autoclass:: lenstools.topology.ConvergenceMap
	:members: fromfilename,gradient,hessian,peakCount,minkowskiFunctionals,moments,powerSpectrum

Shear maps
==========

.. automodule:: lenstools.shear
	:members: load_fits_default

.. autoclass:: lenstools.shear.ShearMap
	:members: fromfilename,decompose,sticks,convergence