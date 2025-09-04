import numpy as np
from pytest import approx

from gnssir.processing import snr2arc


def test_snr2arc_basic(fake_sine_wave_snr_data):
    """Test basic functionality of snr2arc with fake sine wave data."""
    # Use fake data from fixture
    snr_data = fake_sine_wave_snr_data

    # Set up parameters
    rhlims = [2, 4]  # Reflector height limits in meters
    gsignal = "L1"
    polydeg = 2

    # Run snr2arc
    rh_arr, snrdt_arr = snr2arc(snr_data, rhlims, gsignal=gsignal, polydeg=polydeg, detrend=True)

    assert len(rh_arr) == 1

    # Basic checks
    assert isinstance(rh_arr, np.ndarray), "rh_arr should be numpy array"
    assert isinstance(snrdt_arr, np.ndarray), "snrdt_arr should be numpy array"

    assert rh_arr.shape[1] == 12, "rh_arr should have 12 columns"
    # Check that reflector height estimates are within expected range
    heights = rh_arr[:, 1].astype(float)
    assert heights[0] == approx(3, abs=0.1)

    assert snrdt_arr.shape[1] == 5, "snrdt_arr should have 5 columns"


def test_snr2arc_minimal_data(minimal_snr_data):
    """Test snr2arc with minimal data (edge case)."""
    snr_data = minimal_snr_data
    rhlims = [1.0, 5.0]

    # Should handle minimal data gracefully
    rh_arr, snrdt_arr = snr2arc(snr_data, rhlims, detrend=False)

    assert len(rh_arr) == 0
    assert len(snrdt_arr) == 0
