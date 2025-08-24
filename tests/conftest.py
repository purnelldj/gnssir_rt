import pickle

from pytest import fixture


def pload(fname):
    with open(fname, "rb") as f:
        return pickle.load(f)


@fixture
def snr_in():
    return pload("tests/data/snr_in.pkl")


@fixture
def snr_elv_interp_out():
    return pload("tests/data/snr_elv_interp_out.pkl")
