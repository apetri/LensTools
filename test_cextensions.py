#!/usr/bin/env python

"""
Standalone test to verify C extensions work correctly with NumPy 2.x and pytest.
This test doesn't require external data files.
"""

import pytest
import numpy as np
import lenstools
from lenstools.extern import _topology, _gadget2, _nbody, _pixelize


def test_import_c_extensions():
    """Test that all core C extensions can be imported."""
    # These should all succeed
    assert _topology is not None
    assert _gadget2 is not None  
    assert _nbody is not None
    assert _pixelize is not None
    print("✓ All core C extensions imported successfully")


def test_topology_gradient():
    """Test _topology.gradient function with NumPy 2.x."""
    # Create test data
    test_data = np.random.rand(50, 50).astype(np.float64)
    
    # Call gradient function
    grad_x, grad_y = _topology.gradient(test_data, None, None)
    
    # Check results
    assert grad_x.shape == test_data.shape
    assert grad_y.shape == test_data.shape
    assert grad_x.dtype == np.float64
    assert grad_y.dtype == np.float64
    print("✓ _topology.gradient working with NumPy", np.__version__)


def test_numpy_compatibility():
    """Test that C extensions work with current NumPy version (1.x or 2.x)."""
    # Check NumPy version and report it
    major_version = int(np.__version__.split('.')[0])
    print(f"✓ Testing with NumPy {np.__version__} (major version {major_version})")
    
    # Both NumPy 1.x and 2.x should work with our compatibility fixes
    assert major_version >= 1, f"NumPy version too old: {np.__version__}"
    
    # Test array operations with C extensions
    test_array = np.ones((20, 20), dtype=np.float64)
    result = _topology.gradient(test_array, None, None)
    
    # Since gradient of constant array should be zero
    assert np.allclose(result[0], 0.0, atol=1e-10)
    assert np.allclose(result[1], 0.0, atol=1e-10)
    print(f"✓ NumPy {major_version}.x compatibility verified")


def test_lenstools_basic_functionality():
    """Test basic LensTools functionality without external data."""
    # Test that we can create basic objects
    from lenstools.image.convergence import ConvergenceMap
    import astropy.units as u
    
    # Create a simple convergence map
    test_data = np.random.rand(32, 32)
    conv_map = ConvergenceMap(test_data, 1.0 * u.deg)
    
    assert conv_map.data.shape == (32, 32)
    assert conv_map.side_angle.value == 1.0
    print("✓ Basic LensTools functionality working")


if __name__ == "__main__":
    # Run tests directly
    test_import_c_extensions()
    test_topology_gradient() 
    test_numpy_compatibility()
    test_lenstools_basic_functionality()
    print(f"\n🎉 All tests passed! C extensions working with NumPy {np.__version__}")