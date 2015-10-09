.. _catalog_to_map::

Construct convergence maps out of shear catalogs
================================================

Constructing a convergence map of a particular subfield from a shear catalog is an operation that occurs frequently in weak lensing analysis; if the field of view is unmasked, the reconstruction operation is equivalent to the calculation of the :math:`E` mode of the shear field :math:`\gamma`. Here is how to to use lenstools to perform this operation: the operations are handled with the :py:class:`~lenstools.catalog.shear.ShearCatalog` class

::

	from lenstools.catalog.shear import ShearCatalog
	import astropy.table as tbl
	import astropy.units as u

Suppose that the shear table is contained in a file called 'WLshear.fits' (which should contain a table with columns 'shear1','shear2' readable with astropy.table.Table.read()), in which each row represents a different galaxy; suppose also that the sky positioning information :math:`(x,y)` is contained in a file 'positions.fits'; we must combine the information in these two tables to build a :math:`(x,y,\gamma)` table

::

	shear_catalog = ShearCatalog.read('WLshear.fits')
	positions_catalog = ShearCatalog.read('positions.fits')
	full_catalog = tbl.hstack((positions_catalog,shear_catalog))

You can tweak the names of the columns that contain the position information (which by default are 'x','y') and the measure units of the positions (which by default are degrees) using the :py:func:`~lenstools.catalog.shear.ShearCatalog.setSpatialInfo` method. At this point you can construct a shear grid out of the catalog, specifying a grid size and a smoothing scale.

::

	shear_map = full_catalog.toMap(map_size=3.5*u.deg,npixel=512,smooth=0.5*u.arcmin)
	convergence_map = shear_map.convergence()

Note that 'shear_map' now is a :py:class:`~lenstools.image.shear.ShearMap` instance and 'convergence_map' is a :py:class:`~lenstools.image.convergence.ConvergenceMap` instance, and as such you can visualize it with the :py:func:`~lenstools.image.convergence.ConvergenceMap.visualize` method

::
	
	convergence_map.visualize(colorbar=True,cbar_label=r"$\kappa$")

.. figure:: ../../../examples/catalog_to_convergence.png
