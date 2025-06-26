# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

LensTools is a scientific Python package for **Weak Gravitational Lensing analysis** used in astrophysics and cosmology research. It provides tools for convergence maps, shear analysis, N-body simulations, ray tracing, and statistical analysis of cosmic structure.

## Development Commands

### Installation and Build
```bash
# Install in development mode with C extensions
python setup.py build_ext -i --gsl=/usr/local --fftw3=/usr/local
python setup.py install --gsl=/usr/local --fftw3=/usr/local

# Install Python dependencies
pip install -r requirements.txt
```

### Testing
```bash
# Run full test suite with coverage
nosetests --with-coverage --cover-package=lenstools --logging-level=INFO

# Run specific test file
nosetests lenstools/tests/test_convergence.py -v
```

### External Dependencies
The package requires these external C libraries:
- **GSL** (GNU Scientific Library) - for mathematical functions
- **FFTW3** - for fast Fourier transforms  
- **NICAEA** (optional) - for cosmological calculations

Configure library paths in `setup.cfg` if not in `/usr/local`.

## Architecture

### Core Package Structure
- **`image/`** - Convergence maps, shear maps, flexion maps, CMB analysis
- **`statistics/`** - Ensemble analysis, constraints, contours, MCMC samplers
- **`simulations/`** - N-body simulation I/O (Gadget2, FastPM), ray tracing engine
- **`catalog/`** - Galaxy catalog handling and mock generation
- **`pipeline/`** - Simulation pipeline management, cluster computing, MPI support
- **`observations/`** - Real survey data handling (CFHT lens survey)
- **`utils/`** - FFT utilities, MPI helpers, algorithms
- **`extern/`** - C extensions for performance-critical computations
- **`tests/`** - Comprehensive test suite using nose framework

### Key Classes and Concepts
- **`ConvergenceMap`** - Primary class for weak lensing convergence field analysis
- **`ShearMap`/`FlexionMap`** - Shear and higher-order lensing field handling
- **`Catalog`** classes - Galaxy catalog management with astropy integration
- **`Ensemble`** - Statistical ensemble analysis framework for parameter estimation
- **`RayTracer`** - Ray tracing through simulated dark matter structures
- **`Design`** - Optimal experimental design for cosmological surveys

### External Integrations
- **N-body simulations:** Gadget2, FastPM format support
- **Cosmological codes:** CAMB, NICAEA parameter calculations  
- **Observational data:** CFHT lens survey, DES-like mock catalogs
- **Parallel computing:** MPI support for cluster environments

### Performance Features
- C extensions in `extern/` for computationally intensive operations (topology analysis, Gadget I/O, pixelization)
- MPI parallelization throughout the codebase
- Efficient FFT implementations via FFTW3
- NumPy/SciPy vectorized operations

## Testing Framework

Uses **nose** testing framework with coverage reporting. Tests are organized by functionality and located in `lenstools/tests/`. The CI system (Travis) runs tests on Python 3.6/3.7 with external library dependencies.