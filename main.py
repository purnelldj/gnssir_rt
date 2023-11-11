"""
written by David Purnell
https://github.com/purnelldj
"""

from importlib import import_module
import argparse
from gnssr.processing import snr2arcs
from gnssr.processing import arcsplot
from gnssr.processing import arcs2splines

parser = argparse.ArgumentParser()
parser.add_argument("station", help="station ID")
parser.add_argument("funcname", help="function name: snr2arcs, arcsplot or arcs2spline")
parser.add_argument("-t", "--true", nargs='+', help="add [multiple] optional params to = True")
parser.add_argument("-f", "--false", nargs='+', help="add [multiple] optional params to = False")
args = parser.parse_args()

station_id = args.station
funcname = args.funcname

print(station_id, funcname)

station_input_file = "site_inputs." + station_id
try:
    tmod = import_module(station_input_file)
except:
    print(f"station input file {station_input_file} does not exist")
pyargs = tmod.__dict__

# setting optional args True/False
if args.true is not None:
    for arg in args.true:
        pyargs[arg] = True
if args.false is not None:
    for arg in args.false:
        pyargs[arg] = False

if funcname == 'snr2arcs':
    snr2arcs(**pyargs)
elif funcname == 'arcsplot':
    arcsplot(**pyargs)
elif funcname == 'arcs2splines':
    arcs2splines(**pyargs)
else:
    print('did not recognise input function')