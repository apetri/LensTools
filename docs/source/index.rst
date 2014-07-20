.. LensTools documentation master file, created by
   sphinx-quickstart on Sun Jun 29 15:31:39 2014.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to LensTools's documentation!
=====================================

This python package collects together a suite of widely used analysis tools in Weak Gravitational Lensing

.. warning::

	This is still in pre-alpha stage, not tested yet on full scale! Use at your own risk!

Summary
-------

This python add-on will handle basically every operation you will need to perform on Weak Lensing survey data; the distribution includes a range of tools in image analysis, statistical processing and numerical theory predictions and supports multiprocessing using the mpi4py module. 

Installation
------------

The easiest way will be to install through pip::

	pip install LensTools

In the meantime you can clone or fork fork the `github repository <https://github.com/apetri/LensTools>`_ to download the source. The core features require numpy and scipy to run, but the test suite requires additionally the matplotlib and astropy packages. If you want to run the calculations in parallel on a computer cluster you will need to install mpi4py (a python wrapper for MPI) and emcee (from which lenstools borrows the MPI Pool)

Test
----

To check that everything works you can run some pre implemented tests. First you need to install `pytest <http://pytest.org/latest/>`_, then you need to go into the Test directory and run::
	
	py.test

The tests will need some data files to run, the tarball is available `here <https://doc-0k-28-docs.googleusercontent.com/docs/securesc/2onj0g7d6rhqc0ucq5nr982bdvss3iqf/ivlpbbeo0m6u8me3s2a0n34d7462dk07/1405893600000/16262236432539457434/16262236432539457434/0Bx6ZTqNbWBx_Mk01WVJzYm5jNVk?e=download&h=16653014193614665626&nonce=33s37c4aqa36c&user=16262236432539457434&hash=t82hd2amumcn6ti16g7tancug56nc3c8>`_. Watch out, the file is approximately 100MB! After downloading the tarball, unpack it in the Test directory and you are good to go! 

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

