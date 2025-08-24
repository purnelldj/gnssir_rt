# GNSS-IR Real-Time Water Level Processing

[![License: BSD 3-Clause](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

A Python-based command line tool for real-time processing of low-cost GNSS-IR (Global Navigation Satellite System Interferometric Reflectometry) water level data.

## About

This software accompanies the research article ["Real-time water levels using GNSS-IR: a potential tool for flood monitoring"](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2023GL105039) by Purnell et al., published in Geophysical Research Letters (2024).

## Installation

### From PyPI (recommended)

```bash
pip install gnssir-rt
```

### From Source

For an editable install:
```bash
git clone https://github.com/purnelldj/gnssir-rt.git
cd gnssir-rt
pip install -e .
```

## Quick Start

Here's a complete example using test data from Saint-Joseph-de-la-Rive:

### Step 1: Download Test Data

```bash
# Download the test dataset
wget https://github.com/purnelldj/gnssir-rt/raw/main/tests/data/sjdlr.zip

# Create data directory and extract to it
mkdir -p gnssir-data
unzip sjdlr.zip -d gnssir-data
```

### Step 2: Process the Data

```bash
gnssir --task arcs2splines --site_dir gnssir-data/sjdlr --config gnssir-data/sjdlr/sjdlr.yaml
```

This will generate a plot showing one day of processed water level data.

## Usage

### General Command Structure

```bash
gnssir --task TASK --site_dir SITE_DIR --config path/to/config.yaml
```

### Configuration Options

You can either change parameters directly in a configuration file or from the command line as shown below.

```bash
gnssir --task arcs2splines --site_dir gnssir-data/sjdlr --config gnssir-data/sjdlr/sjdlr.yaml \
       --azilims '[190, 230]'
```

### Help and Documentation

```bash
gnssir --help  # Show all available parameters
```

## Data

### Research Dataset

The complete SNR dataset used in the accompanying research paper is available on Zenodo: [10.5281/zenodo.10114719](https://doi.org/10.5281/zenodo.10114719)

To download the dataset:
```bash
pip install zenodo_get
zenodo_get 10.5281/zenodo.10114719
```

**Note:** The download is approximately 1GB.

### SNR Data Format

The software expects SNR data in a specific format similar to the [gnssrefl format](https://gnssrefl.readthedocs.io/en/latest/pages/file_structure.html#the-snr-data-format), with the following modifications:

| Column | Description |
|--------|-------------|
| 1 | Satellite number |
| 2 | Elevation angle (degrees) |
| 3 | Azimuth angle (degrees) |
| 4 | [GPS time](https://docs.astropy.org/en/stable/api/astropy.time.TimeGPS.html) (instead of seconds of day) |
| 5 | L1 SNR values |

## Contributing

We welcome contributions from the community! Whether you're fixing bugs, adding features, improving documentation, or sharing new site configurations, your help is appreciated.

### How to Contribute

1. **Fork the repository** on GitHub
2. **Clone your fork** locally:
   ```bash
   git clone https://github.com/yourusername/gnssir-rt.git
   cd gnssir-rt
   ```
3. **Create a new branch** for your feature:
   ```bash
   git checkout -b feature/your-feature-name
   ```
4. **Make your changes** and test them thoroughly
5. **Commit your changes** with clear, descriptive messages:
   ```bash
   git commit -m "Add feature: description of your changes"
   ```
6. **Push to your fork** and submit a pull request

### Development Setup

For development work:
```bash
pip install -e ".[dev]"  # Install with development dependencies
python -m pytest         # Run tests
```

## License

This project is licensed under the BSD 3-Clause License - see the [LICENSE](LICENSE) file for details.

## Citation

If you use this software in your research, please cite:

```
Purnell, D., et al. (2024). Real-time water levels using GNSS-IR: a potential tool
for flood monitoring. Geophysical Research Letters.
https://doi.org/10.1029/2023GL105039
```

## Acknowledgments

- The script [make_gpt.py](./gnssir/make_gpt.py) was adapted from the [gnssrefl](https://github.com/kristinemlarson/gnssrefl) repository
- This work builds upon the broader GNSS-IR research community's contributions

## Support

- **Issues**: Report bugs or request features via [GitHub Issues](https://github.com/purnelldj/gnssir-rt/issues)
- **Discussions**: Join community discussions in [GitHub Discussions](https://github.com/purnelldj/gnssir-rt/discussions)
- **Email**: For direct inquiries, contact the maintainers

---

**Keywords:** GNSS-IR, water level monitoring, flood detection, interferometric reflectometry, remote sensing
