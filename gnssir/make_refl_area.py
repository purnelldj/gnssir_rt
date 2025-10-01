"""
This code has been adapted from an earlier version of 'gnssrefl'
https://github.com/kristinemlarson/gnssrefl/

Functions for calculating and visualizing GNSS-IR reflection areas.
"""

import numpy as np
import simplekml
from pathlib import Path
from typing import List, Tuple, Optional

from gnssir.tropd import corr_rh_facs


def get_fresnel_dimensions(rh: float, elv: float, gsignal: str = "L1") -> List[float]:
    """
    Calculate Fresnel zone dimensions for GNSS-IR.

    Parameters
    ----------
    rh : float
        Reflector height in meters
    elv : float
        Elevation angle in degrees
    gsignal : str, optional
        GNSS signal type, by default "L1"

    Returns
    -------
    List[float]
        [R, a, b] where R is radial distance to center, a and b are semi-axes

    Raises
    ------
    ValueError
        If gsignal is not supported
    """
    if gsignal == "L1":
        lfreq = 1575.42e6
    else:
        raise ValueError("Only L1 signal is currently supported")

    lcar = 299792458 / lfreq  # wavelength
    n = 1  # first fresnel zone - don't change
    d = n * lcar / 2

    sin_elv = np.sin(np.radians(elv))
    tan_elv = np.tan(np.radians(elv))

    R = rh / tan_elv + (d / sin_elv) / tan_elv  # radial dist to center
    b = np.sqrt(2 * d * rh / sin_elv + (d / sin_elv) ** 2)
    a = b / sin_elv

    return [R, a, b]


def calculate_circle_lat_lon(
    lat: float, lon: float, radius: float, azimuth_limits: List[float], num_points: int = 10
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Calculate latitude and longitude points for a circle arc on Earth's surface.

    Parameters
    ----------
    lat : float
        Center latitude in degrees
    lon : float
        Center longitude in degrees
    radius : float
        Radius in meters
    azimuth_limits : List[float]
        [min_azimuth, max_azimuth] in degrees
    num_points : int, optional
        Number of points to generate, by default 10

    Returns
    -------
    Tuple[np.ndarray, np.ndarray]
        (latitudes, longitudes) arrays of circle points

    Raises
    ------
    ValueError
        If azimuth limits are reversed
    """
    earth_radius = 6.371e6  # radius of earth in meters

    if azimuth_limits[1] < azimuth_limits[0]:
        raise ValueError("Upper and lower azimuth limits are reversed")

    azi = np.linspace(azimuth_limits[0], azimuth_limits[1], num_points + 1)
    azi_rad = np.radians(azi)
    ang_dist = radius / earth_radius  # angular distance

    sin_lat = np.sin(np.radians(lat))
    cos_lat = np.cos(np.radians(lat))
    lon_rad = np.radians(lon)

    lat_circ = np.arcsin(sin_lat * np.cos(ang_dist) + cos_lat * np.sin(ang_dist) * np.cos(azi_rad))
    lon_circ = lon_rad + np.arctan2(
        np.sin(azi_rad) * np.sin(ang_dist) * cos_lat,
        np.cos(ang_dist) - sin_lat * np.sin(lat_circ),
    )

    lat_circ = np.degrees(lat_circ)
    lon_circ = np.degrees(lon_circ)

    return lat_circ, lon_circ


def calculate_reflection_area(
    lla: List[float],
    rh_limits: List[float],
    azimuth_limits: List[float],
    elevation_limits: List[float],
    gsignal: str = "L1",
    tropd_adj: bool = True,
    full_fresnel: bool = False,
    station_name: Optional[str] = None,
    output_dir: str = ".",
) -> dict:
    """
    Calculate and generate KML file for GNSS-IR reflection area.

    Parameters
    ----------
    lla : List[float]
        [latitude, longitude, altitude] of the station
    rh_limits : List[float]
        [min_height, max_height] reflector height limits in meters
    azimuth_limits : List[float]
        [min_azimuth, max_azimuth] in degrees
    elevation_limits : List[float]
        [min_elevation, max_elevation] in degrees
    gsignal : str, optional
        GNSS signal type, by default "L1"
    tropd_adj : bool, optional
        Apply tropospheric delay adjustment, by default True
    full_fresnel : bool, optional
        Use full Fresnel zone, by default False
    station_name : Optional[str], optional
        Name for the reflection area, by default None
    output_dir : str, optional
        Output directory for KML file, by default "."

    Returns
    -------
    dict
        Dictionary containing polygon coordinates and metadata
    """
    lat_lon = [lla[0], lla[1]]

    if tropd_adj:
        # Default atmospheric parameters (not critical for area calculation)
        pant = 1000  # pressure in hPa
        tant = 10  # temperature in °C
        eant = 5  # water vapor pressure in hPa

        rhfac_outer = corr_rh_facs(elevation_limits[0], elevation_limits[0], pant, tant, eant)
        rhfac_inner = corr_rh_facs(elevation_limits[1], elevation_limits[1], pant, tant, eant)

        ffz_outer = get_fresnel_dimensions(
            rh_limits[1] + rh_limits[1] * rhfac_outer, elevation_limits[0], gsignal=gsignal
        )
        ffz_inner = get_fresnel_dimensions(
            rh_limits[0] + rh_limits[0] * rhfac_inner, elevation_limits[1], gsignal=gsignal
        )
    else:
        ffz_outer = get_fresnel_dimensions(rh_limits[1], elevation_limits[0], gsignal=gsignal)
        ffz_inner = get_fresnel_dimensions(rh_limits[0], elevation_limits[1], gsignal=gsignal)

    rad_outer = ffz_outer[0]
    rad_inner = ffz_inner[0]

    if full_fresnel:
        rad_outer = rad_outer + ffz_outer[1]
        rad_inner = rad_inner - ffz_inner[1]

    # Calculate circle coordinates
    lat_outer, lon_outer = calculate_circle_lat_lon(
        lat_lon[0], lat_lon[1], rad_outer, azimuth_limits, num_points=10
    )
    lat_inner, lon_inner = calculate_circle_lat_lon(
        lat_lon[0], lat_lon[1], rad_inner, azimuth_limits, num_points=5
    )

    # Create polygon by going backwards along inner circle, then forwards along outer circle
    lat_poly = np.append(np.flip(lat_inner), lat_outer)
    lat_poly = np.append(lat_poly, lat_inner[-1])
    lon_poly = np.append(np.flip(lon_inner), lon_outer)
    lon_poly = np.append(lon_poly, lon_inner[-1])

    # Create coordinate pairs for different formats
    lon_lat_pairs = [[lon_poly[i], lat_poly[i]] for i in range(len(lat_poly))]  # for GEE

    # Create WKT polygon string
    snap_poly = "POLYGON(("
    for i in range(len(lon_poly)):
        if i > 0:
            snap_poly += ", "
        coord = f"{lon_poly[i]:.10f} {lat_poly[i]:.10f}"
        snap_poly += coord
    snap_poly += "))"
    print(f"WKT Polygon:\n{snap_poly}")

    # Create KML file
    kml = simplekml.Kml()
    refl_name = station_name if station_name else "refl_area"

    pol = kml.newpolygon(name=refl_name, outerboundaryis=lon_lat_pairs)
    pol.style.polystyle.color = simplekml.Color.changealphaint(150, simplekml.Color.pink)

    # Ensure output directory exists
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    kml_path = f"{output_dir}/{refl_name}.kml"
    kml.save(kml_path)
    print(f"\nKML file saved to: {kml_path}")

    return {
        "gee_coordinates": lon_lat_pairs,
        "wkt_polygon": snap_poly,
        "kml_path": kml_path,
    }
