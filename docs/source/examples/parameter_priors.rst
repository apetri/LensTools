A quick and dirty way to incorporate parameter priors
=====================================================

Suppose you carried on your weak lensing analysis all the way to the parameter constraints, and you were able to estimate your parameter covariance matrix :math:`\Sigma_{lens}` (either from simulated or real data). Now suppose you are interested in understanding how these constraints change when you add prior information from say CMB observations from Planck. These prior results will become available through their own parameter covariance matrix :math:`\Sigma_{CMB}`, which may, or may not, have the same dimensions and parametrization as :math:`\Sigma_{lens}`. Applying the prior to the parameters considered in the weak lensing analysis and fixing all the others is equivalent to take the appropriate parameter slice of :math:`\Sigma_{CMB}^{-1}` and adding the Fisher matrices

.. math:: \Sigma_{lens+CMB} = (\Sigma_{lens}^{-1}+\Sigma_{CMB}^{-1})^{-1}

This can be readily done with the functionality embedded in the :py:class:`~lenstools.statistics.ensemble.SquareMatrix` class, with the following code

::

	from lenstools.statistics.ensemble import SquareMatrix

	#Read in parameter covariances
	lens_pcov = SquareMatrix.read("lenscov.pkl")
	cmb_cov = SquareMatrix.read("cmbcov.pkl")

	#Parametrization
	parameters = ["Om","w","sigma8"]

	#Add the Fisher matrices
	fisher_lens_cmb = lens_pcov.invert()[parameters] + cmb_cov.invert()[parameters]

	#pcov_lens_cmb is the parameter covariance subject to the prior
	pcov_lens_cmb = fisher_lens_cmb.invert()
