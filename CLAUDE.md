# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

LensTools is a scientific Python package for **Weak Gravitational Lensing analysis** used in astrophysics and cosmology research. It provides tools for convergence maps, shear analysis, N-body simulations, ray tracing, and statistical analysis of cosmic structure.

## Development Commands

### Installation and Build
```bash
# Modern installation (preferred method)
pip install -e .

# Or with dependencies from requirements file
pip install -r requirements.txt
pip install -e .

# Legacy installation with C extensions (if GSL/FFTW3 available)
python setup.py build_ext -i --gsl=/usr/local --fftw3=/usr/local  
python setup.py install --gsl=/usr/local --fftw3=/usr/local
```

### Testing
```bash
# Modern testing with pytest (recommended)
pytest lenstools/tests/ -v

# Legacy testing with nose (requires nose package, incompatible with Python 3.13+)
nosetests --with-coverage --cover-package=lenstools --logging-level=INFO

# Test basic functionality and C extensions
python -c "import lenstools; from lenstools.extern import _topology; print('All working!')"
```

### Modernization Status ✅ COMPLETE

The package has been successfully modernized with the following improvements:

- **Modern packaging**: Added `pyproject.toml` with modern build system
- **NumPy 2.x compatibility**: ✅ COMPLETE - All C extensions now work with NumPy 2.0+
- **Updated dependencies**: Modern versions of scipy, matplotlib, astropy, emcee
- **C extensions**: ✅ WORKING - All C extensions compile and run correctly
- **Cross-platform Python support**: Python 3.8+ compatibility
- **Legacy compatibility**: Maintains compatibility with older Python versions

**C Extensions Status**: 
- ✅ All C extensions (_topology, _gadget2, _nbody, _pixelize, _design, _nicaea) are enabled and working
- ✅ NumPy 2.x C API compatibility issues resolved using compatibility macros
- ✅ All PyArray_DATA, PyArray_DIM calls updated for NumPy 2.x

**Remaining Minor Issues**:
- Some deprecated NumPy aliases (np.float, np.complex) need updating throughout codebase
- pkg_resources warnings can be addressed by migrating to importlib.resources

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