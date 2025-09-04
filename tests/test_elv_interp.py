import pickle
from pathlib import Path

from pytest import approx

from gnssir.elv_interp import elv_interp_array


def load_test_data(filename):
    """Load test data from pickle file."""
    test_data_dir = Path(__file__).parent / "data"
    with open(test_data_dir / filename, "rb") as f:
        return pickle.load(f)


def test_elv_interp_array():
    """Test elevation interpolation with real data."""
    # Load test data
    snr_in = load_test_data("snr_in.pkl")
    snr_elv_interp_out = load_test_data("snr_elv_interp_out.pkl")

    # Test the function
    snr_out = elv_interp_array(snr_in, kspac=1800)

    # Assertions
    assert snr_in.shape[1] == snr_out.shape[1], "mismatching shapes from elv interp"
    assert approx(snr_elv_interp_out) == snr_out, "mismatching values from elv interp"
