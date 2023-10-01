# -*- coding: utf-8 -*-
"""
Authors: Erik Hedenström and Anton Petrén, 2022
Part of the Bachelor's thesis "Using Jupiter's Moon Io as a Plasma Probe"

This file contains functions for obtaining the coordinates of a celestial body
at a given time. It is necessary to download appropriate Kernels to be able to
run this code.
"""

import spiceypy
from datetime import datetime, timedelta
import numpy as np

jupiterRadius = 71492  # jupiter radius in km


def getPosition(object, t):
    """Returns the position of the given celestial body in jovicentric
    cartesian coordinates"""
    if object == Europa:
        cart = Spice(Europa, t)
    if object == Io:
        cart = Spice(Io, t)
    if object == Jupiter:
        cart

    return cart


def cartesianToPolar(cart):
    """Return [r, theta, phi] in jupiter radii and degrees (theta=0 at equator)"""
    
    x = cart[0]
    y = cart[1]
    z = cart[2]
    r = np.linalg.norm(cart)
    theta = np.arccos(z/r)*180/np.pi
    phi = np.arctan(y/x)*180/np.pi
    if x < 0:
        phi += 180
    if phi < 0:
        phi += 360
    theta -= 90
    theta = theta*(-1)
    return [r/jupiterRadius,theta,phi]


def datetimeToET(t):
    """Convert datetime object to ET time"""
    
    utctim = t.strftime("%Y %B %d %H:%M:%S")
    return spiceypy.str2et(utctim)  # Convert utctime to ET.


def Spice(object, t):
    """Uses spice to get a set of coordinates (vector) in jovicentric
    cartesian coordinates for the celestial body (string). t=datetime object"""
    
    METAKR = 'Kernels.tm'  # name of file with kernels to use
    spiceypy.furnsh(METAKR)
    et = datetimeToET(t)
    [pos, ltime] = spiceypy.spkpos(object, et, 'IAU_JUPITER', 'LT+S', 'Jupiter')
    return pos


def localTime(object, t):
    """Return local datetime at earthtime t for object (string)"""
    
    METAKR = 'Kernels.tm'  # name of file with kernels to use
    spiceypy.furnsh(METAKR)
    et = datetimeToET(t)
    [pos, ltime] = spiceypy.spkpos('EARTH', et, 'J2000', 'LT+S', object, )
    tloc = t - timedelta(seconds=int(ltime))
    return tloc


def main():
    # used only for testing
    d = datetime(year=2015, month=6, day=11, hour=0, minute=0, second=0)
    print("Time on Earth:", d.strftime("%Y %B %d %H:%M:%S"))
    print("Io position:", Spice("Io", d))
    print("Io position, polar:", cartesianToPolar(Spice("Io", d)))
    print("Europa position:", Spice("Europa", d))
    print("Sun position:", Spice("Sun", d))  # light travel time?
    print("Jupiter position:", Spice("Jupiter", d))
    print("Local time on Jupiter:", localTime("Jupiter", d).strftime("%Y %B %d %H:%M:%S"))


if __name__ == "__main__":
    main()
