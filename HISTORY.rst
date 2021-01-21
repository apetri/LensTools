.. :changelog:

1.2
+++

- Support for weak lensing flexion analysis (courtesy of Brij Patel, brp53@drexel.edu) 

1.0
+++

- Plot Gaussian prediction for peak histograms (Bond et. al. 1987)
- Make reduced shear maps and catalogs in lenstools.raytracing
- Better handling of random seeds in raytracing
- Measure equilateral and folded bispectrum of scalar images (kappa,cmb)
- Introduced CMBTemperatureMap class
- CMB lensing potential reconstructions with quadratic TT estimator (wrapped around quicklens)
- Handle neutrino parameters in pipeline book keeping
- Introduced Gadget2SnapshotNu class for N-body outputs with neutrino effects

0.9
+++

- RayTracer can now perform line of sight integrations to approximate kappa (Born, lens-lens and geodesic perturbation) and omega (lens-lens), without full raytracing
- cutPlanes can now handle snapshots in FastPM format
- If snapshots are provided in light cone projections, LensTools can do plane cutting, raytracing and LOS integrations in a single script
- Fit multiple observations with the same FisherAnalysis  

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


