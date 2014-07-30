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

This python add-on will handle basically every operation you will need to perform on Weak Lensing survey data; the distribution includes a range of tools in image analysis, statistical processing and numerical theory predictions and supports multiprocessing using the `mpi4py <http://mpi4py.scipy.org>`_ module. 

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

The core features require `numpy <http://www.numpy.org>`_ and `scipy <http://www.scipy.org>`_ to run, and the Test suite requires additionally the `matplotlib <http://matplotlib.org>`_ and `astropy <http://www.astropy.org>`_ packages. If you want to run the calculations in parallel on a computer cluster you will need to install `mpi4py <http://mpi4py.scipy.org>`_ (a python wrapper for the MPI library) and `emcee <http://dan.iel.fm/emcee/current/>`_ (from which LensTools borrows the MPI Pool utility)

Test
----

To check that everything works before installing you can run the pre implemented test suite that comes with the source code. First you will need to install `pytest <http://pytest.org/latest/>`_, then you need to download some data files (mainly FITS images) that the test suite depends on (you can find these data files `here <https://doc-0k-28-docs.googleusercontent.com/docs/securesc/2onj0g7d6rhqc0ucq5nr982bdvss3iqf/ivlpbbeo0m6u8me3s2a0n34d7462dk07/1405893600000/16262236432539457434/16262236432539457434/0Bx6ZTqNbWBx_Mk01WVJzYm5jNVk?e=download&h=16653014193614665626&nonce=33s37c4aqa36c&user=16262236432539457434&hash=t82hd2amumcn6ti16g7tancug56nc3c8>`_, the file is approximately 100MB). Once you downloaded it into the Test directory unpack it::

   tar -xvf data.tar.gz

Then you're good to go! Just run::
	
	py.test

Each test, if successful, will produce some PNG plots.

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

Gallery
-------

.. toctree::
   :maxdepth: 2

   gallery

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

