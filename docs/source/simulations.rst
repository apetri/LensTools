Weak lensing simulations
************************

At the current stage of Weak gravitational lensing research, large numerical simulations are required for analyzing observations; the LensTools python package provides an API to interact with some already existing simulated datasets (mainly convergence and shear maps for different cosmological models), such as 

1. The IGS1 simulations: this simulated dataset contains 1000 realizations of single redshift convergence and shear maps for six different cosmological parammeter combinations (a fiducial model and some variations). The fiducial model is based on 45 independent N-body simulations and the variations are based on 5 independent N-body simulations (where :math:`N=512^3`)

2. The CFHTemu1 simulations: this simulated dataset contains 1000 realizations of convergence maps with the source redshift distribution of the CFHTLens survey; the simulated foregrounds are available for 91 different cosmological parameter variations of the triplet :math:`(\Omega_m,w,\sigma_8)`

While LensTools provides the API to interact with a local copy of these simulated datasets, it does not provide the datasets themselves; to obtain the IGS1 dataset please email `Jan M. Kratochvil <jan.m.kratochvil@gmail.com>`_, to obtain the CFHTemu1 dataset please email `Andrea Petri <apetri@phys.columbia.edu>`_