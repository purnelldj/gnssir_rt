# gnssir_rt
This software is to go with "Real-time water levels using GNSS-IR: a potential tool for flood monitoring" by Purnell et al. (under review for AGU GRL)
All code written by David Purnell except for gnssr/make_gpt.py (written by Kristine Larson)

## Dependencies
numpy, astropy, matplotlib, scipy

## How to use the code
SNR data can be coverted to 'arcs' and then to a water level spline.
Easiest to use from the command line as follows
```
python main.py [station] [funcname]
```
where `[station]` corresponds to a file: `site_inputs/[station].py`
and funcname is one of 'snr2arcs', 'arcsplot' or 'arcs2splines'

## SNR data format
SNR data has been provided to go with the paper, it can be found at [link]
The SNR data format is the same as: https://gnssrefl.readthedocs.io/en/latest/pages/file_structure.html#the-snr-data-format
but with differences on fourth and fifth+ columns:
Instead of seconds of day in the fourth column it is GPS time: https://docs.astropy.org/en/stable/api/astropy.time.TimeGPS.html
The fifth column is L1 SNR (there are only five columns)
