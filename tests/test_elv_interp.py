import pickle

from pytest import approx

from gnssir_rt.elv_interp import elv_interp_array


def pload(fname):
    with open(fname, "rb") as f:
        return pickle.load(f)


# taking 1hr from one antenna
def test_elv_interp_array():
    snr_in = pload("tests/testdata/snr_elv_interp_in.pkl")
    snr_out_saved = pload("tests/testdata/snr_elv_interp_out.pkl")
    snr_out = elv_interp_array(snr_in, kspac=1800)
    assert snr_in.shape[1] == snr_out.shape[1]
    assert approx(snr_out_saved) == snr_out
