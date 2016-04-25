.. :changelog:

0.7-dev
+++++++

- ShearCatalog allows to re-bin galaxies according to redshift
- Introduced a SquareMaxtrix class (inherits from Ensemble) for square matrix operations, with column name support
- Protect nodes in a SimulationBatch to call forbidden methods

0.6
+++

- Database class for local/remote SQL databases
- Gadget2SnapshotPipe class allows to combine Gadget2 and lenstools.planes via named pipes
- Improvements in pipeline deployment 

0.5
+++

- Weak lensing pipeline functionalities
- Measure 2pcf of convergence maps with hankel transforms
- Moved Limber module under lenstools.simulations
- Dedicated module for Fast Fourier Transforms
- Ensemble is now a sub-class of pandas.DataFrame


0.4.8.4
+++++++

- New operations in Ensemble: bootstrap, readall, shuffle
- Ensemble mean is automatically recomputed when re--assigning data

0.4.8
+++++

Beta release 


