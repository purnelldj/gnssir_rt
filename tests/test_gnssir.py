import pickle
from importlib import import_module
from os import listdir
from shutil import rmtree

from gnssir_rt.processing import arcs2splines, snr2arcs


def pload(fname):
    with open(fname, "rb") as f:
        return pickle.load(f)


# 1. test some arc file outputs
def test_arc_out():
    station_input_file = "tests.sjdlr_t1input"
    tmod = import_module(station_input_file)
    pyargs = tmod.__dict__
    arcdirbm = pyargs["arcdir"]
    arcdirtest = arcdirbm + "_temp"
    pyargs["arcdir"] = arcdirtest
    snr2arcs(**pyargs)
    for aid in pyargs["antennaids"]:
        antdir = arcdirtest + "/" + aid
        testfile = listdir(antdir)
        for testf in testfile:
            fulltestfile = antdir + "/" + testf
            fullbmfile = arcdirbm + "/" + aid + "/" + testf
            testarr = pload(fulltestfile)
            bmarr = pload(fullbmfile)
            assert (testarr == bmarr).all()
    # now remove arcdirtest
    rmtree(arcdirtest)


# 2. test one day of output spline
def test_spline():
    station_input_file = "tests.sjdlr_t2input"
    tmod = import_module(station_input_file)
    pyargs = tmod.__dict__
    pyargs["plotfig"] = False
    invout_test = arcs2splines(**pyargs)
    invoutbmfile = pyargs["stationdir"] + "/invout.pkl"
    invout_bm = pload(invoutbmfile)
    for key in invout_bm:
        assert (invout_bm[key] == invout_test[key]).all()
    # below was for creating the bm
    # with open(invoutbmfile, 'wb') as f:
    #    pickle.dump(invout, f)


# 3. test tropd input file
# skipping for now - not my code
