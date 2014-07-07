LensTools tests
================

To test the package on your machine you will need to install the python module test package (pytest); just type in a shell

    pip install pytest

After that you will need a test map to run the tests on: you will need to put in this directory a file called _conv.fit_ which should be a FITS image which has a keyword _ANGLE_ in the header (one of the convergence maps we have been using for example); if you want the test to succeed also for the shear utilities you should add two shear maps too named _shear1.fit_,_shear2.fit_ After that just run

    py.test

This will test if the code runs and if everything succeeds a bunch of png figures will appear. Good luck!

P.S. Depending on how you installed the package, you may or may not need to run

    python setup.py build_ext -i

to compile the C backends!!
