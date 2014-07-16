.. LensTools documentation master file, created by
   sphinx-quickstart on Sun Jun 29 15:31:39 2014.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to LensTools's documentation!
=====================================

This python package collects together a suite of widely used analysis tools in Weak Gravitational Lensing

Summary
-------

This python add-on will handle basically every operation you will need to perform on Weak Lensing survey data; the distribution includes a range of tools in image analysis, statistical processing and numerical theory predictions and supports multiprocessing using the mpi4py module. 

Installation
------------

The easiest way is to install through pip::

	pip install LensTools

Test
----

To check that everything works you can run some pre implemented tests. First you need to install `pytest <http://pytest.org/latest/>`_, then you need to go into the Test directory and run::
	
	py.test

The tests will need some data files to run: you can download these from ...

Issues
------

The code is maintained and developed on `github <https://github.com/apetri/LensTools>`_, pull requests are welcome!

License
-------

None so far

API
---

.. toctree::
   :maxdepth: 2

   code

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

