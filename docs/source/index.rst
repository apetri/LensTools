.. LensTools documentation master file, created by
   sphinx-quickstart on Sun Jun 29 15:31:39 2014.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to LensTools!
=====================

This python package collects together a suite of widely used analysis tools in Weak Gravitational Lensing

Summary
-------

This python add-on will handle basically every operation you will need to perform on Weak Lensing survey data; the distribution includes a range of tools in image analysis, statistical processing and numerical theory predictions and supports multiprocessing using the mpi4py_ module. This package includes many useful features, including:

- Measure units handling through astropy
- Complete flexibility and easy customization of input/output formats
- Efficient measurements of power spectrum, PDF, Minkowski functionals and peak counts of convergence maps
- Survey masks
- Efficient :math:`E` and :math:`B` mode decompositions of shear maps
- Artificial noise generation engines
- All functionality of pandas_ DataFrame readily available (and improved!)
- Easy to compute parameter statistical inferences
- Easy input/output from :math:`N`-body simulation snapshots in the Gadget2 binary format
- Interfaces with existing simulation sets
- Ray Tracing simulations
- CPU vectorization of expensive computations via numpy
- Easy multiprocessing and cluster deployment via the mpi4py_ module
- *Future prospect*: taking advantage of numpy offload capabilities to Intel Xeon Phi coprocessors to boost performance (planned) 

Acknowledgement
---------------

LensTools is open source and released under the MIT license. If you make use of LensTools in your work, please cite it!

::
   
	@ARTICLE{2016A&C....17...73P,
   		author = {{Petri}, A.},
    		title = "{Mocking the weak lensing universe: The LensTools Python computing package}",
  		journal = {Astronomy and Computing},
		archivePrefix = "arXiv",
   		eprint = {1606.01903},
 		keywords = {Weak Gravitational Lensing, Simulations},
     		year = 2016,
    		month = oct,
   		volume = 17,
    		pages = {73-79},
      		doi = {10.1016/j.ascom.2016.06.001},
   		adsurl = {http://adsabs.harvard.edu/abs/2016A%26C....17...73P},
  		adsnote = {Provided by the SAO/NASA Astrophysics Data System}
	}

The above code paper is published on the `Elsevier Astronomy and Computing journal <http://www.journals.elsevier.com/astronomy-and-computing/>`_. 

Installation
------------


The easiest way is to install through pip ::

	pip install lenstools

The LensTools installer will look for the optional dependencies (GSL, FFTW3, NICAEA) in the standard location /usr/local. If you want to specity a different location for these dependencies you need to specify it with the --install-option switch. For example running ::

   pip install lenstools --install-option="--gsl=/my/path/to/gsl" 

the installer will look for GSL under /my/path/to/gsl. If you do this don't forget to update the LD_LIBRARY_PATH environment variable ::

   export LD_LIBRARY_PATH=/my/path/to/gsl:$LD_LIBRARY_PATH

to ensure a correct dynamic linking for the external libraries. 
An alternative is to install from source by cloning or forking the `github repository <https://github.com/apetri/LensTools>`_ to download the source and build it manually. First clone the repository (the original one, or your fork)::
   
   git clone https://github.com/apetri/LensTools

Then, inside the LensTools directory build the source::

   python setup.py build

If you want to test the build before installing to your system, look to the instructions in Test. Once you are satisfied install the package (might require root privileges depending on the install location)::

   python setup.py install

Dependencies
------------

.. _numpy: http://www.numpy.org
.. _scipy: http://www.scipy.org
.. _astropy: http://www.astropy.org
.. _pandas: http://pandas.pydata.org
.. _emcee: http://dan.iel.fm/emcee/current/
.. _matplotlib: http://matplotlib.org
.. _mpi4py: http://mpi4py.scipy.org
.. _GSL: http://www.gnu.org/software/gsl/    
.. _pytest: http://pytest.org/latest/
.. _here: http://www.plus-six.com/andrea/lenstools/data.tar.gz
.. _github: https://github.com/apetri/LensTools
.. _NICAEA: http://www.cosmostat.org/software/nicaea/
.. _fftw3: http://www.fftw.org

The core features require the standard numpy_, scipy_ , and additionally  astropy_ (mainly for the cosmology and measure units support) and emcee_ (from which LensTools borrows the MPI Pool utility), and the Test suite requires additionally the matplotlib_ package. matpoltlib should eventually be installed if you want to use the plotting engines of LensTools. If you want to run the calculations in parallel on a computer cluster you will need to install mpi4py_ (a python wrapper for the MPI library). Installation of all these packages is advised (if you run astrophysical data analyses you should use them anyway). One of the lenstools features, namely the :py:class:`~lenstools.simulations.Design` class, requires that you have a working version of GSL_ to link to; if you don't have one, just hit *enter* during the installation process and the package will work correctly without this additional feature. The installation if the NICAEA_ bindings additionally requires a working installation of the fftw3_ library. 

Test
----

To check that everything works before installing you can run the pre implemented test suite that comes with the source code. First you will need to install pytest_, then you need to download some data files (mainly FITS images) that the test suite depends on. You need to set the environment variable LENSTOOLS_DATA to the path where you want your data to be downloaded (for a manual download the data file can be found here_, it is almost 250MB). After that, in a python shell, type

::

   import lenstools
   lenstools.dataExtern()

That should make sure the data directory is downloaded and available and should return the full path of the data directory. After that operation has completed create a Test directory and run the tests:

::
   
   mkdir Test
   cd Test
   py.test --pyargs lenstools.tests

Each test, if successful, will produce some output (PNG plots, text files, etc...)


Weak Lensing simulations
------------------------

.. toctree::
   :maxdepth: 2

   pipeline
   raytracing
   simulations

Command line scripts
--------------------

.. toctree::
   :maxdepth: 2

   scripts

Gallery
-------

.. toctree::
   :maxdepth: 2

   gallery

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

Issues
------

The code is maintained and developed on github_, pull requests are welcome!

License
-------

Copyright 2014 Andrea Petri and contributors; lenstools is free software made available under the MIT License. For details see the LICENSE file


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

