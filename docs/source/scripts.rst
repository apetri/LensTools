LensTools command line scripts
******************************

General purpose scripts
=======================

nbodyheader
-----------

Displays the header information of a Nbody simulation snapshot. Usage:

::
	
	nbodyheader <file_1> <file_2> ...



gadgetstrip
-----------

Strips the sometimes unnecessary velocity information from a Gadget2 snapshot, roughly reducing its disk usage by a factor of 2. Usage:

::
	
	gadgetstrip <file_1> <file_2> ...

The script can also read from standard input: to strip all the files in a directory:

::
	
	ls | gadgetstrip


npyinfo
-------

lenstools.confidencecontour
---------------------------

lenstools.view
--------------


LensTools pipeline scripts
==========================

lenstools.submission
--------------------

lenstools.cutplanes
-------------------

lenstools.cutplanes-mpi
-----------------------

lenstools.raytracing
--------------------

lenstools.raytracing-mpi
------------------------

lenstools.execute-mpi
---------------------