.. LensTools documentation master file, created by
   sphinx-quickstart on Sun Jun 29 15:31:39 2014.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to LensTools!
=====================

This python package collects together a suite of widely used analysis tools in Weak Gravitational Lensing

.. warning::

	This is still in pre-alpha stage, not tested yet on full scale! Use at your own risk!

Summary
-------

This python add-on will handle basically every operation you will need to perform on Weak Lensing survey data; the distribution includes a range of tools in image analysis, statistical processing and numerical theory predictions and supports multiprocessing using the `mpi4py <http://mpi4py.scipy.org>`_ module. This package includes many useful features, including:

- Measure units handling through astropy
- Complete flexibility and easy customization of input/output formats
- Efficient measurements of power spectrum, PDF, Minkowski functionals and peak counts of convergence maps
- Survey masks
- Efficient :math:`E` and :math:`B` mode decompositions of shear maps
- Artificial noise generation engines
- Easy to compute parameter statistical inferences
- Easy input/output from :math:`N`-body simulation snapshots in the Gadget2 binary format
- Interfaces with existing simulation sets
- Ray Tracing simulations (in development...)
- CPU vectorization of expensive computations via numpy
- Easy multiprocessing and cluster deployment via the mpi4py module 

Installation
------------

The easiest way WILL be to install through pip (coming really soon, not available yet)::

	pip install LensTools

In the meantime you can clone or fork the `github repository <https://github.com/apetri/LensTools>`_ to download the source and build it manually. 
First clone the repository (the original one, or your fork)::
   
   git clone https://github.com/apetri/LensTools

Then, inside the LensTools directory build the source::

   python setup.py build

If you want to test the build before installing to your system, look to the instructions in Test. Once you are satisfied install the package (might require root privileges depending on the install location)::

   python setup.py install

Dependencies
------------

The core features require `numpy <http://www.numpy.org>`_, `scipy <http://www.scipy.org>`_ and `astropy <http://www.astropy.org>`_ (mainly for the cosmology and measure units support)  to run, and the Test suite requires additionally the `matplotlib <http://matplotlib.org>`_ package. matpoltlib should eventually be installed if you want to use the plotting engines of LensTools. If you want to run the calculations in parallel on a computer cluster you will need to install `mpi4py <http://mpi4py.scipy.org>`_ (a python wrapper for the MPI library) and `emcee <http://dan.iel.fm/emcee/current/>`_ (from which LensTools borrows the MPI Pool utility). Installation of all these packages is advised (if you run astrophysical data analyses you should use them anyway)

Test
----

To check that everything works before installing you can run the pre implemented test suite that comes with the source code. First you will need to install `pytest <http://pytest.org/latest/>`_, then you need to download some data files (mainly FITS images) that the test suite depends on (you can find these data files `here <http://danishlaundromat.com/apetri/data.tar.gz>`_, the file is approximately 250MB). Once you downloaded it into the Test directory unpack it::

   tar -xvf data.tar.gz

Then you're good to go! Just run::
	
	py.test

Each test, if successful, will produce some PNG plots.

Issues
------

The code is maintained and developed on `github <https://github.com/apetri/LensTools>`_, pull requests are welcome!

License
-------

Copyright 2014 Andrea Petri and contributors; lenstools is free software made available under the MIT License. For details see the LICENSE file



Gallery
-------

.. toctree::
   :maxdepth: 2

   gallery

Weak lensing simulations
------------------------

.. toctree::
   :maxdepth: 2

   simulations
   raytracing

3D visualization with Mayavi
----------------------------

.. toctree::
   :maxdepth: 2

   visualization

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

