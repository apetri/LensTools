LensTools tests
================

To test the package on your machine you will need to install the python module test package (pytest); just type in a shell

    pip install pytest

After that you will need a test map to run the tests on: you will need to put in this directory a file called _map.fit_ which should be a FITS image which has a keyword _ANGLE_ in the header (one of the convergence maps we have been using for example). After that just run

    py.test

This will test if the code runs and if everything succeeds a bunch of png figures will appear. Good luck!