import datetime
import os

# find your local directory
gnssir_rt_dir = os.getcwd()
gnssir_rt_dir = gnssir_rt_dir.replace("\\", "/")  # for windows

# general parameters
stationdir = gnssir_rt_dir + "/localproc/sjdlr"
sdt = datetime.datetime(2021, 11, 25, 0)  # zenodo data range: 11, 10, 22 to 11, 27, 18
edt = datetime.datetime(2021, 11, 26, 0)  #
antennaids = ["ACM0", "ACM1", "ACM2", "ACM3"]
hgts = [0.2, 0.3, 0, 0.1]
lla = [47.4488045, -70.365557, -20]
gptfile = stationdir + "/tropd_input.txt"
snrdir = stationdir + "/snr"
arcdir = stationdir + "/arcs"

# snr2arc parameters
elvinterp = True
azilims = [190, 250]  # 190, 250
elvlims = [5, 20]  # 5, 20
rhlims = [1.5, 9]
tempres = 5
arclim = 1 * 60 * 60
iterdt_arcs = 1 * 60 * 60
snrfilelen = 1 * 60 * 60

# arc2spline parameters
kdt = 2 * 60 * 60
iterdt_spline = 1 * 60 * 60
prek = 2
postk = 2
fixed_std = 1
# latency = 60*60
pktnlim = 3

outdir = "temp_outputs"
