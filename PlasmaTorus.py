# -*- coding: utf-8 -*-
"""
Authors: Erik Hedenström and Anton Petrén, 2022
Part of the Bachelor's thesis "Using Jupiter's Moon Io as a Plasma Probe"

This file is used for calculating the plasma density at a given location based
on the model of the plasma and the parameters provided in the Test file.
"""

import numpy as np
import math
import scipy.integrate as integrate
from mpmath import csc, cot
# R = number of jupiter radii


def plasmaDensity(coordinates, flat = False, integrate = False,
                 N=[1710,2180,2160,1601],
                 H=[0.1,0.6,1.0],
                 W=[0.20,0.08,0.32,1.88],
                 C=[5.23,5.60,5.89,5.53],
                 C2min=5.6, C2max=5.76,
                 circle=True, phi0=0, theta0f=0):
    """Calculates the plasma density at a given set of coordinates"""

    eqDev, theta0 = thetaToEqDev(coordinates=coordinates, flat=flat, integrate=integrate, phi0=phi0)

    R = coordinates[0]
    theta = coordinates[1]
    phi = coordinates[2]
    sunPhi = coordinates[3]

    if integrate:
        R = fieldLineR(R, theta0, theta)
    density = mathFunction(eqDev, R, phi, sunPhi, N, H, W, C, C2min, C2max, circle, phi0)
    return density, eqDev


def fieldLineR(R, theta0, theta):
    theta = math.radians(theta)
    theta0 = math.radians(theta0)
    if np.cos(theta-theta0)==0:
        return R*1000  # = infinity
    return R/np.cos(theta-theta0)**2


def thetaToEqDev(coordinates, phi0, flat = False, integrate = False):
    """Calculates the distance between the given coordinate and the centrifugal equator"""
    
    R = coordinates[0]
    theta = coordinates[1]
    phi = coordinates[2]

    if not flat:
        theta0 = centrifugalEquator(R, phi, phi0)
    else:
        theta0 = flatCentrifugalEquator(R, phi, phi0)

    if integrate: # follow magnetic field lines
        eqDev = integratingFunction(R, theta, theta0)
    else: # vertical torus
        eqDev = R*(math.sin(math.radians(theta))-math.sin(math.radians(theta0)))

    return eqDev, theta0


def centrifugalEquator(R, phi, phi0):
    """Returns the angle between the equatorial plane and the centrifugal equator
    for a curved model of the centrifugal equator"""
    
    #calculates in degrees but calculates sin and cos in radians
    phi = phi - phi0

    a = 1.66        #degrees
    b = 0.131
    c = 1.62
    d = 7.76        #degrees
    e = 249         #degrees
    theta0 = (a*math.tanh(math.radians(b*R-c))+d)*math.sin(math.radians(phi-e)) #returned as lattitude from jupiters equator
    return theta0


def flatCentrifugalEquator(R, phi, phi0):
    """Returns the angle between the equatorial plane and the centrifugal equator,
    if the centrifugal equator plane is assumed to be flat"""
    
    phi = phi - phi0
    e = 249  #degrees
    theta0 = 6.4*math.sin(math.radians(phi-e)) #returned as latitude from jupiters equator
    return theta0


def integratingFunction(R, theta, theta0):
    """Returns distance from the centrifugal equator to a coordinate, integrated
    along the magnetic field lines"""
    
    theta = math.radians(90-theta)
    theta0 = math.radians(90-theta0)
    theta00 = math.radians(theta0)
    op = integrate.quad(toIntegrate, theta, theta0, args = (theta, theta00, R))
    return op[0]


def toIntegrate(x, a, c, r):
    """Dipole magnetic field"""
    
    op = r * math.sqrt( (4*cot(x+c)**2 + 1)*csc(x+c)**4 )
    return op


def magneticEquator(phi, theta0f):
    """Returns the tilt of the magnetic equator for given phi"""
    
    e = 249  # degrees
    theta0 = theta0f + 6.4*math.sin(math.radians(phi - e))  # returned as latitude from jupiters equator
    return theta0


def isDusk(phi, sunPhi, phi0=0):
    """Return True if phi is on dusk side, False otherwise"""
    
    phiDusk = sunPhi + 90 + phi0  # phi0 = deviation from real dusk direction
    v = (phi - phiDusk) * np.pi / 180  # IO relative dusk side
    i = np.cos(v)
    return i > 0


def calcC(C, C2min, C2max, phi, sunPhi, circle, phi0):
    """Return new C with respect to dusk-dawn effects"""
    
    dC2 = C2max - C2min
    
    # calculate dusk direction
    phiDusk = sunPhi + 90 + phi0  # phi0 = deviation from real dusk direction
    
    # calculate where Io is relative the dusk direction
    v = (phi-phiDusk)*np.pi/180

    if circle:
        C2 = C2min + 0.5*dC2*(-np.cos(v) + 1)
    else:
        ec = (C2max-C2min)/(C2max+C2min)
        p = C2max*(1-ec)
        C2 = p / (1 - -ec*np.cos(v))

    # move the whole torus together with the ribbon
    dC = C2 - C[1]
    C1 = C[0] + dC
    C3 = C[2] + dC
    C4 = C[3] + dC
    C = [C1, C2, C3, C4]
    return C


def mathFunction(eqDev, R, phi, sunPhi, N, H, W, C, C2min, C2max, circle=True, phi0=0):
    """Uses the known mathematical functions to calculate plasma density
    according to the model of the plasma torus"""
    
    N1 = N[0]
    N2 = N[1]
    N3 = N[2]
    N4 = N[3]
    H1 = H[0]
    H2 = H[1]
    H3 = H[2]
    W1 = W[0]
    W2 = W[1]
    W3 = W[2]
    W4 = W[3]
    C20 = C[1]

    if C2max - C2min > 0:
        # Implementing the dawn-dusk effect into C
        C = calcC(C, C2min, C2max, phi, sunPhi, circle, phi0)

    # peak location (Jupiter radii):
    C1 = C[0]  # cold torus
    C2 = C[1]  # ribbon
    C3 = C[2]  # warm torus
    C4 = C[3]  # extended torus

    dC = C2 - C20
    if R < 6.1 + dC:
        density = N1 * np.exp(-(R - C1) ** 2 / W1 ** 2) * np.exp(-eqDev ** 2 / H1 ** 2) + \
                  N2 * np.exp(-(R - C2) ** 2 / W2 ** 2) * np.exp(-eqDev ** 2 / H2 ** 2) + \
                  N3 * np.exp(-(R - C3) ** 2 / W3 ** 2) * np.exp(-eqDev ** 2 / H3 ** 2)
    else:
        density = N4 * np.exp(-(R - C4) ** 2 / W4 ** 2) * np.exp(-eqDev ** 2 / H3 ** 2)

    return density


def main():
    # used only for testing
    print(isDusk(234.09783817008685, 85.89092893018517))  # brightest data on dusk
    coordinates = [5.901375204971329, 0.003571640236231133, 234.09783817008685, 85.89092893018517]
    plasmaDensity(coordinates, C2min=5.56, C2max=5.84, circle=True)

    print(isDusk(139.66106665365308, 314.19285577638857))  # dimmest data on dawn
    coordinates = [5.8809787905452735, 0.005312128221902412, 139.66106665365308, 314.19285577638857]
    plasmaDensity(coordinates, C2min=5.56, C2max=5.84, circle=True)


if __name__ == "__main__":
    main()
