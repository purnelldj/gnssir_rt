import numpy as np
from pytest import fixture


@fixture
def fake_sine_wave_snr_data():
    """Create fake SNR data based on a sine wave pattern for testing."""
    # Create fake SNR data with sine wave pattern
    # Columns: [sat_prn, elevation, azimuth, time, snr]

    n_points = 100
    sat_prn = 5  # GPS satellite
    time_start = 1000.0  # GPS time
    time_end = 2000.0
    elv_start = 5.0  # degrees
    elv_end = 25.0  # degrees
    azimuth = 45.0  # degrees (constant)

    # Create time and elevation arrays
    times = np.linspace(time_start, time_end, n_points)
    elevations = np.linspace(elv_start, elv_end, n_points)
    azimuths = np.full(n_points, azimuth)
    sat_prns = np.full(n_points, sat_prn)

    # Create sine wave SNR pattern based on elevation
    # Simulate reflector height oscillations
    reflector_height = 3.0  # meters
    wavelength = 0.19  # L1 wavelength in meters
    frequency = 2 * reflector_height / wavelength

    # Convert elevation to sine of elevation for SNR calculation
    sin_elv = np.sin(elevations * np.pi / 180)

    # Create sine wave pattern in SNR (simulating multipath interference)
    base_snr = 45.0  # dB-Hz
    amplitude = 5.0  # dB-Hz
    phase_factor = 2 * np.pi * frequency
    snr_values = base_snr + amplitude * np.sin(phase_factor * sin_elv)

    # Add some noise
    np.random.seed(42)  # For reproducible tests
    noise = np.random.normal(0, 0.5, n_points)
    snr_values += noise

    # Stack into required format: [sat_prn, elevation, azimuth, time, snr]
    snr_data = np.column_stack([sat_prns, elevations, azimuths, times, snr_values])

    return snr_data


@fixture
def minimal_snr_data():
    """Create minimal SNR data for edge case testing."""
    # Very simple data with just a few points
    snr_data = np.array(
        [
            [5, 10.0, 45.0, 1000.0, 40.0],  # sat, elv, azi, time, snr
            [5, 15.0, 45.0, 1100.0, 42.0],
            [5, 20.0, 45.0, 1200.0, 41.0],
        ]
    )
    return snr_data
