from datetime import datetime
from typing import Dict, Any
import os
from dataclasses import dataclass
from typing import List

from jsonargparse import ArgumentParser

from gnssir.processing import arcs2splines, arcsplot, snr2arcs


@dataclass
class Config:
    """Configuration for GNSS-IR processing tasks."""

    # required parameters
    task: str  # Processing task: snr2arcs, arcsplot, or arcs2splines
    sdt: str  # Start datetime (YY-MM-DD HH:MM)
    edt: str  # End datetime (YY-MM-DD HH:MM)
    site_dir: str  # Parent directory containing site data in subdirectories 'snr' and 'arcs'
    antennaids: List[str]  # List of antenna IDs (e.g., ["ACM0", "ACM1", "ACM2", "ACM3"])
    hgts: List[float]  # Antenna heights [m] (e.g., [0.2, 0.3, 0, 0.1])
    lla: List[float]  # Site coordinates [lat, lon, height] (e.g., [47.4488045, -70.365557, -20])
    azilims: List[float]  # Azimuth limits [min, max] degrees (e.g., [190, 250])
    elvlims: List[float]  # Elevation limits [min, max] degrees (e.g., [5, 20])
    rhlims: List[float]  # Reflector height limits [min, max] meters (e.g., [1.5, 9])

    # Directories and files
    snr_dir: str = None  # SNR data directory (default: site_dir/snr)
    arc_dir: str = None  # Arc data directory (default: site_dir/arcs)
    gptfile: str = None  # GPT file path (default: site_dir/tropd_input.txt)
    outdir: str = ".gnssir_tmp"  # Output directory

    # SNR to arcs parameters
    elvinterp: bool = True  # Enable elevation interpolation
    tempres: int = 5  # Temporal resolution (minutes)
    arclim: int = 3600  # Arc time limit (seconds)
    iterdt_arcs: int = 3600  # Arc iteration time step (seconds)
    snrfilelen: int = 3600  # SNR file length (seconds)

    # Arcs to splines parameters
    prek: int = 2  # Pre-knot parameter
    postk: int = 2  # Post-knot parameter
    iterdt_spline: int = 3600  # Spline iteration time step (seconds)
    pktnlim: int = 3  # Peak point limit
    kdt: int = 7200  # Knot spacing time (seconds)
    fixed_std: float = 1.0  # Fixed standard deviation for spline fitting


def main(config: Config) -> None:
    """
    Main entry point for GNSS-IR processing tasks.

    This function handles different processing tasks based on the configuration:
    - snr2arcs: Convert SNR data to arcs
    - arcsplot: Plot the arcs
    - arcs2splines: Convert arcs to splines

    Parameters
    ----------
    config : Config
        Configuration object containing processing parameters

    Raises
    ------
    ValueError
        If the specified task is not recognized
    """

    pyargs = load_cfg(config)

    if config.task == "snr2arcs":
        snr2arcs(**pyargs)
    elif config.task == "arcsplot":
        arcsplot(**pyargs)
    elif config.task == "arcs2splines":
        arcs2splines(**pyargs)
    else:
        print(f"input function '{config.task}' not recognized ")
        print("e.g., set --task=snr2arcs from command line")


def load_cfg(config: Config) -> Dict[str, Any]:
    """
    Load and process configuration parameters.

    Converts date strings in the configuration to datetime objects and
    maps parameter names to match the processing functions.

    Parameters
    ----------
    config : Config
        Configuration object

    Returns
    -------
    Dict[str, Any]
        Processed configuration dictionary with datetime objects
    """
    pyargs = config.__dict__.copy()

    # Convert datetime strings
    for dt in ["sdt", "edt"]:
        pyargs[dt] = datetime.strptime(pyargs[dt], "%y-%m-%d %H:%M")
    # make directories and files relative to site_dir if not explicitly provided
    if pyargs["snr_dir"] is None:
        pyargs["snr_dir"] = os.path.join(pyargs["site_dir"], "snr")
    if pyargs["arc_dir"] is None:
        pyargs["arc_dir"] = os.path.join(pyargs["site_dir"], "arcs")
    if pyargs["gptfile"] is None:
        pyargs["gptfile"] = os.path.join(pyargs["site_dir"], "tropd_input.txt")

    for arg in pyargs:
        print(f"{arg}: {pyargs[arg]}")

    return pyargs


def cli() -> None:
    """CLI entry point using jsonargparse."""
    parser = ArgumentParser(
        description="GNSS-IR Real-Time Water Level Processing",
        default_config_files=[],  # Don't auto-load config files
    )
    parser.add_class_arguments(Config)

    # Add config file argument
    parser.add_argument("--config", action="config")

    args = parser.parse_args()
    config = parser.instantiate_classes(args)

    main(config)


if __name__ == "__main__":
    cli()
