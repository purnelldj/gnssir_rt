import numpy as np
import pytest
from pytest import approx

from gnssir.helper import (
    glonasswlen,
    cubspl_nans,
    residuals_cubspl_spectral,
)


def test_glonasswlen():
    """Test GLONASS wavelength calculation."""
    # Test with valid PRN for L1 signal
    prn = 101
    wavelength_l1 = glonasswlen(prn, "L1")
    assert isinstance(wavelength_l1, float)
    assert wavelength_l1 == approx(299792458 / (1602e06 + 0.5625e06))
    # Test with valid PRN for L2 signal
    wavelength_l2 = glonasswlen(prn, "L2")
    assert isinstance(wavelength_l2, float)
    assert wavelength_l2 == approx(299792458.0 / (1246e06 + 0.4375e06))


def test_glonasswlen_invalid():
    """Test GLONASS wavelength calculation with invalid input."""
    # Test with invalid PRN (should raise exception)
    with pytest.raises(Exception):
        glonasswlen(130, "L1")  # PRN > 124 should fail

    # Test with invalid signal (should return NaN)
    result = glonasswlen(101, "L5")
    assert np.isnan(result)


def test_cubspl_nans():
    """Test cubic spline interpolation with NaN handling."""
    # Create simple test data
    knots = np.array([0, 1, 2, 3, 4])
    kval = np.array([1.0, 2.0, 3.0, 4.0, 5.0])  # Linear data
    tplot = np.array([0.5, 1.5, 2.5, 3.5])

    result = cubspl_nans(tplot, knots, kval)
    assert len(result) == len(tplot)
    # For linear data, interpolation should give values close to expected
    expected = np.array([1.5, 2.5, 3.5, 4.5])
    assert approx(result, abs=0.1) == expected


def test_cubspl_nans_with_nans():
    """Test cubic spline interpolation with NaN values in input."""
    # Create test data with NaN
    knots = np.array([0, 1, 2, 3, 4])
    kval = np.array([1.0, np.nan, 3.0, 4.0, 5.0])  # NaN in middle
    tplot = np.array([0.5, 2.5, 3.5])

    result = cubspl_nans(tplot, knots, kval)
    assert np.isnan(result[0])
    expected = np.array([3.5, 4.5])
    assert approx(result[1:], abs=0.1) == expected


def test_residuals_cubspl_spectral():
    """Test spectral residuals calculation."""
    # Create simple test data
    kval = np.array([2.0, 2.5, 3.0, 3.5, 4.0])
    knots = np.array([100, 200, 300, 400, 500])  # Time points

    # Create rh_arr with required columns [time, height, sat, tane_dedt, ...]
    rh_arr = np.array(
        [
            [150, 2.2, 1, 0.1, 10, 20, 45],
            [250, 2.8, 2, 0.15, 15, 25, 90],
            [350, 3.2, 3, 0.12, 12, 22, 135],
        ],
        dtype=object,
    )

    residuals = residuals_cubspl_spectral(kval, knots, rh_arr)
    assert len(residuals) == len(rh_arr)  # Should have one residual per data point
    assert residuals == approx(np.array([0, 0, 0]), abs=0.1)
