"""
This code has been adapted from an earlier version of 'gnssrefl'
https://github.com/kristinemlarson/gnssrefl/
"""

import numpy as np
import simplekml

from gnssir.tropd import corr_rh_facs


def getFresDims(rh, elv, gsignal="L1"):
    if gsignal == "L1":
        lfreq = 1575.42e6
    else:
        raise Exception("only works for L1 right now")
    lcar = 299792458 / lfreq
    n = 1  # first fresnel zone - don't change
    d = n * lcar / 2
    sinElv = np.sin(elv / 180 * np.pi)
    tanElv = np.tan(elv / 180 * np.pi)
    R = rh / tanElv + (d / sinElv) / tanElv  # radial dist to center
    b = np.sqrt(2 * d * rh / sinElv + (d / sinElv) ** 2)
    a = b / sinElv
    ffz = [R, a, b]
    return ffz


def circleLatLon(lat, lon, radius, azilims, NumPoints=10, **kwargs):
    R = 6.371e6  # radius of earth in meters
    if azilims[1] < azilims[0]:
        raise Exception("upper and lower azi lims reversed - exit")
    azi = np.linspace(azilims[0], azilims[1], NumPoints + 1)
    aziRad = [az / 180 * np.pi for az in azi]
    angDist = radius / R  # distance is always radius out from station
    sinlat = np.sin(lat / 180 * np.pi)
    coslat = np.cos(lat / 180 * np.pi)
    lonRad = lon / 180 * np.pi
    latCirc = np.arcsin(sinlat * np.cos(angDist) + coslat * np.sin(angDist) * np.cos(aziRad))
    lonCirc = lonRad + np.arctan2(
        np.sin(aziRad) * np.sin(angDist) * coslat,
        np.cos(angDist) - sinlat * np.sin(latCirc),
    )
    latCirc = latCirc / np.pi * 180
    lonCirc = lonCirc / np.pi * 180
    # plt.subplots(figsize=(5, 5))
    # plt.plot(lonCirc, latCirc)
    # tmp = ''
    # if 'tmp' in kwargs:
    #    tmp = kwargs.get('tmp')
    # plt.savefig(tmp + 'circtest.png', format='png')
    return latCirc, lonCirc


def reflarea(lla, rhlims, azilims, elvlims, gsignal="L1", tropdAdj=True, fullFres=False, **kwargs):
    latLon = [lla[0], lla[1]]
    if tropdAdj:
        pant = 1000  # random average values, not really important
        tant = 10
        eant = 5
        rhfacOuter = corr_rh_facs(elvlims[0], elvlims[0], pant, tant, eant)
        rhfacInner = corr_rh_facs(elvlims[1], elvlims[1], pant, tant, eant)
        ffzOuter = getFresDims(rhlims[1] + rhlims[1] * rhfacOuter, elvlims[0], gsignal=gsignal)
        ffzInner = getFresDims(rhlims[0] + rhlims[0] * rhfacInner, elvlims[1], gsignal=gsignal)
    else:
        ffzOuter = getFresDims(rhlims[1], elvlims[0], gsignal=gsignal)
        ffzInner = getFresDims(rhlims[0], elvlims[1], gsignal=gsignal)
    radOuter = ffzOuter[0]
    radInner = ffzInner[0]
    if fullFres:
        radOuter = radOuter + ffzOuter[1]
        radInner = radInner - ffzInner[1]
    latOuter, lonOuter = circleLatLon(latLon[0], latLon[1], radOuter, azilims, NumPoints=10)
    latInner, lonInner = circleLatLon(latLon[0], latLon[1], radInner, azilims, NumPoints=5)
    # assuming that looking southward, go clockwise
    # go backwards along inner circle
    # then forwards along outter circle
    latPoly = np.append(np.flip(latInner), latOuter)
    latPoly = np.append(latPoly, latInner[-1])
    lonPoly = np.append(np.flip(lonInner), lonOuter)
    lonPoly = np.append(lonPoly, lonInner[-1])
    # lonLatPairs = [(lonPoly[i], latPoly[i]) for i in range(len(latPoly))]
    lonLatPairs = [[lonPoly[i], latPoly[i]] for i in range(len(latPoly))]  # for ee
    lonLatSNAP = [(lonPoly[i], latPoly[i]) for i in range(len(latPoly))]
    # snapPoly = 'POLYGON((%3.10f %3.10f, %3.10f %3.10f, %3.10f %3.10f, %3.10f %3.10f, %3.10f %3.10f))' %(lonInner[-1], latInner[-1],
    # lonInner[0], latInner[0], lonOuter[0], latOuter[0], lonOuter[-1], latOuter[-1], lonInner[-1], latInner[-1])
    print("for GEE:")
    print(lonLatPairs)
    print("for SNAP:")
    print(lonLatSNAP)
    snapPoly = "POLYGON(("
    for i in range(len(lonPoly)):
        if i > 0:
            snapPoly = snapPoly + ", "
        tcoor = "%3.10f %3.10f" % (lonPoly[i], latPoly[i])
        snapPoly = snapPoly + tcoor
    snapPoly = snapPoly + "))"
    print(snapPoly)
    kml = simplekml.Kml()
    reflName = "reflArea"
    if "stationName" in kwargs:
        stationName = kwargs.get("stationName")
        reflName = stationName
    pol = kml.newpolygon(name=reflName, outerboundaryis=lonLatPairs)  # lon, lat, optional height
    pol.style.polystyle.color = simplekml.Color.changealphaint(150, simplekml.Color.pink)  # max alpha is 255
    if "tmp" in kwargs:
        tmp = kwargs.get("tmp")
    kml.save(tmp + "reflArea.kml")
    return
