# gnssir_rt
This software is to go with "Real-time water levels using GNSS-IR: a potential tool for flood monitoring" by Purnell et al. (Geophysical Research Letters)
All code written by David Purnell except for gnssr/make_gpt.py (written by Kristine Larson)

## Installation
Requires python >= 3.9, git and pip.

### Step 1: clone repository
```
git clone https://github.com/purnelldj/gnssir_rt
cd gnssir_rt

```
### Step 2: [optional] create and activate venv
If you skip this step then the dependencies will be installed on your global python environment. To create a venv called '.venv' using mac OS:
```
python -m venv .venv
source .venv/bin/activate

```
### Step 3: install using pip
```
python -m pip install .
```


## How to use the code
SNR data can be coverted to 'arcs' and then to a water level spline.
Easiest to use from the command line as follows
```
python main.py [station] [funcname]
```
where `[station]` corresponds to a file: `site_inputs/[station].py`
and funcname is one of 'snr2arcs', 'arcsplot' or 'arcs2splines'

## SNR data format
SNR data to go with the paper can be found [here](https://doi.org/10.5281/zenodo.10114719). The SNR data format is similar to [this format](https://gnssrefl.readthedocs.io/en/latest/pages/file_structure.html#the-snr-data-format), but with differences on fourth and fifth columns:
* instead of seconds of day in the fourth column it is [GPS time](https://docs.astropy.org/en/stable/api/astropy.time.TimeGPS.html)
* the fifth column is L1 SNR (there are only five columns)
