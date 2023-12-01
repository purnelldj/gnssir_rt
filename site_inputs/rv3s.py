import datetime
import os

# find your local directory
gnssir_rt_dir = os.getcwd()
gnssir_rt_dir = gnssir_rt_dir.replace("\\", "/")  # for windows

# general parameters
stationdir = gnssir_rt_dir + "/localproc/rv3s"
sdt = datetime.datetime(2020, 9, 12, 0)  # zenodo date range 9, 9, 18 to 10, 9, 13
edt = datetime.datetime(2020, 9, 13, 0)  #
antennaids = ["a", "b", "c", "d"]
hgts = [0.2, 0.3, 0, 0.1]
lla = [46.340526, -72.539128, -22.4]
gptfile = stationdir + "/tropd_input.txt"
snrdir = stationdir + "/snr"
arcdir = stationdir + "/arcs"

# snr2arc parameters
elvinterp = True
azilims = [80, 220]  # 80, 220
elvlims = [5, 30]  # 10, 20
rhlims = [3, 7]
tempres = 5
arclim = 1 * 60 * 60
iterdt_arcs = 1 * 60 * 60
snrfilelen = 1 * 60 * 60

# arc2spline parameters
kdt = 1 * 60 * 60
iterdt_spline = 1 * 60 * 60
prek = 2
postk = 2
fixed_std = 0.3
# latency = 60*60
pktnlim = 3

outdir = "temp_outputs"
