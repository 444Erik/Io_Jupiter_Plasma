# -*- coding: utf-8 -*-
"""
Author: Erik Hedenström (2023)

This file is used for plotting the plasma torus. It might take a few minutes
to produce a plot (depending on resolution and computer performance).
"""

import ObjectPosition
import PlasmaTorus
import numpy as np
import multiprocessing as mp
from datetime import datetime, timedelta
import os
import re
from mayavi import mlab
from astropy.coordinates import solar_system_ephemeris
from tvtk.api import tvtk
solar_system_ephemeris.set("jpl")

jupiterRadius = 71492  # jupiter radius in km

# Number of simultaneous multithreading processes
NUMBEROFPROCESSES = 8

# Plasma torus model parameters
N = [1710, 2180, 2160, 1601]
H = [0.1, 0.6, 1.0]
W = [0.20, 0.08, 0.32, 1.88]
C = [5.23, 5.60, 5.89, 5.53]  # Default values for torus with constant radius
circle = True  # True for circular torus, False for elliptic torus
fieldLines = False  # True for a torus that follows the magnetic field lines, False for a vertical torus
flat = False  # True for a flat centrifugal equator, False for a curved centrifugal equator
theta0f = 0  # Extra magnetic field tilt (degrees). NOT IMPLEMENTED. It is
# currently assumed that the magnetic equator coincides with the centrifugal equator

#Dawn-dusk effects
C2min = 5.56  # Ribbon minimum distance from center of Jupiter. Overrides C[1]
C2max = 5.84  # Ribbon maximum distance from center of Jupiter. Overrides C[1]
phi0 = 0  # Deviation from true dawn/dusk-direction (degrees)
# It is assumed that C2max is on the dusk side and that C2min is on the dawn
# side. The value provided in C[1] above should be the average ribbon distance,
# smaller than C2max and larger than C2min. All other C values will be updated
# such that the whole torus moves with the position of the ribbon.
# If C2min=C2max, dawn-dusk effects will be ignored and the C[1] value above
# will be used for the ribbon.

# SPICE options
startTime = datetime(year=1980, month=8, day=23, hour=19, minute=15, second=0)
endTime = datetime(year=2022, month=8, day=24, hour=1, minute=0, second=0)
moon = "Io"

# Plotting options
DATE = datetime(year=1998, month=8, day=23, hour=21, minute=1, second=30) - timedelta(seconds=2015.815)
PHISUN = ObjectPosition.cartesianToPolar(ObjectPosition.Spice("Sun", DATE))[1]
RESOLUTION = 100
DATAPOINTS = RESOLUTION**3+1
RMAX = 8
ZMAX = 5

# Choose what to plot in the "runProcessesModified" function

# Save density data
docname = 'plasmatorus2.npy'
if os.path.exists(docname):
    os.remove(docname) # deleting any previous version of this document
# Note that this can be imported and plotted in the "main" function


def coordinates(n):
    """Return list of n cartesian coordinates in km"""
    
    # km, cartesian
    R = jupiterRadius
    dn = int(DATAPOINTS**(1/3))
    X = np.linspace(-RMAX*R, RMAX*R, dn)
    Y = np.linspace(-RMAX*R, RMAX*R, dn)
    Z = np.linspace(-ZMAX*R, ZMAX*R, dn)
    coords = []
    for x in X:
        for y in Y:
            for z in Z:
                coords.append([x,y,z])
    return coords
    

def runProcessesModified(obj, flat, integrate, dusk, time, N, H, W, C, C2min, C2max, circle, phi0, theta0f):
    """This main function serves only as a function to start processes and graph the data they return, all calculation are done in other pocesses"""
    
    with mp.Manager() as manager:
        processes = []
        lock = mp.Lock()
        
        data = coordinates(DATAPOINTS)
        results = manager.list([]) #results = coordinates and density

        for i in range(NUMBEROFPROCESSES):
            p = mp.Process(target = testChain, args = (obj, data, lock, i, NUMBEROFPROCESSES, results, flat, integrate, N, H, W, C, C2min, C2max, circle, phi0, theta0f))
            processes.append(p)

        for p in processes:
            p.start()

        for p in processes:
            p.join()
            
        print("Number of data points:", len(results))
        #mathematica(results)
        fig = mlab.figure(size=(1200, 900), bgcolor=(0,0,0))
        plotTorus(fig, results)
        plotJupiter(fig)
        mlab.show()
        #return results


def testChain(obj, data, lock, i, n, results, flat, integrate, N, H, W, C, C2min, C2max, circle, phi0, theta0f):
    """Allows us to use a single process to run multiple tests, lowering the minimum number of processes"""
    
    for x in range(i, len(data), n):
        result = singleTestModified(obj, data[x], flat, integrate, N, H, W, C, C2min, C2max, circle, phi0, theta0f)
        lock.acquire()
        results.append((data[x], result)) # append coordinates and density
        lock.release()


def singleTestModified(obj, coordinate, flat, integrate, N, H, W, C, C2min, C2max, circle, phi0, theta0f):
    """Runs the code we wish to use to test a single data point and returns the value we wish to display"""
    
    posPolar = ObjectPosition.cartesianToPolar(coordinate)
    arg = [posPolar[0], posPolar[1], posPolar[2], PHISUN]
    density, eqDev = PlasmaTorus.plasmaDensity(coordinates=arg, flat=flat, integrate=integrate, N=N, H=H, W=W, C=C, C2min=C2min, C2max=C2max, circle=circle, phi0=phi0, theta0f=theta0f)

    return density


def plotTorus(fig, results):
    """Plot the plasma torus in Python"""
    
    print("Start python plot function")
    resultsStr = []
    for r in results:
        # append all results to a string to create a dictionary
        resultsStr.append((str([int(r[0][0]),int(r[0][1]),int(r[0][2])]), r[1]))
    resultsDict = dict(resultsStr)

    # km, cartesian
    R = jupiterRadius
    dn = int(DATAPOINTS**(1/3))
    X = np.linspace(-RMAX*R, RMAX*R, dn)
    Y = np.linspace(-RMAX*R, RMAX*R, dn)
    Z = np.linspace(-ZMAX*R, ZMAX*R, dn)

    n = RESOLUTION
    D = np.zeros((n,n,n))

    for xi in range(len(X)):
        for yi in range(len(Y)):
            for zi in range(len(Z)):
                x = int(X[xi])
                y = int(Y[yi])
                z = int(Z[zi])
                D[xi,yi,zi] = resultsDict['[{}, {}, {}]'.format(x,y,z)]  # add density of the correct coordinate to the correct position in D
                
    # Save density data to a file
    with open(docname, 'wb') as f:
        np.save(f, D)
    
    mlab.pipeline.volume( mlab.pipeline.scalar_field(D, colormap='autumn'))  # , vmin=0.3, vmax=0.6 
    return


def plotJupiter(fig):
    """Add Jupiter to the plot"""
    
    max_distance = RMAX*jupiterRadius
    grid_size = max_distance*2  # km
    n = RESOLUTION
    pixel_size = grid_size/n  # km
    print("Pixel size (km):", pixel_size)
    scale = jupiterRadius/pixel_size
    x = (n+1)/2
    y = (n+1)/2
    z = (n+1)/2
    
    # load and map the texture
    image_file = 'jupiter.jpg'
    img = tvtk.JPEGReader()
    img.file_name = image_file
    texture = tvtk.Texture(input_connection=img.output_port, interpolate=1)
    # (interpolate for a less raster appearance when zoomed in)
    
    # use a TexturedSphereSource
    R = scale
    Nrad = 180
    
    # create the sphere source with a given radius and angular resolution
    sphere = tvtk.TexturedSphereSource(radius=R, theta_resolution=Nrad, phi_resolution=Nrad)
    
    # assemble rest of the pipeline, assign texture    
    sphere_mapper = tvtk.PolyDataMapper(input_connection=sphere.output_port)
    sphere_actor = tvtk.Actor(mapper=sphere_mapper, texture=texture, position=(x,y,z))
    fig.scene.add_actor(sphere_actor)
    return


def plotSphere(fig):
    """Add a sphere to the plot"""

    max_distance = RMAX*jupiterRadius
    grid_size = max_distance*2  # km
    n = RESOLUTION
    pixel_size = grid_size/n  # km
    print("Pixel size:", pixel_size)
    scale = jupiterRadius/pixel_size*2
    x = (n+1)/2
    y = (n+1)/2
    z = (n+1)/2
    mlab.points3d(x, y, z, scale_factor=scale, resolution=50, color=(0.75,0.1,0))
    fig
    return


def mathematica(results):
    """Write result to a file in mathematica syntax"""
    
    # Save density data in mathematica format
    docname = 'plasmatorus.m'
    
    print("Start mathematica function")
    resultsStr = []
    for r in results:
        # append all results to a string
        resultsStr.append((str([int(r[0][0]),int(r[0][1]),int(r[0][2])]), r[1]))
    resultsDict = dict(resultsStr)

    X = []
    Y = []
    Z = []

    for i in range(len(results)):
        c = results[i][0]
        x = int(c[0])
        y = int(c[1])
        z = int(c[2])
        if x not in X:
            X.append(x)
        if y not in Y:
            Y.append(y)
        if z not in Z:
            Z.append(z)

    X.sort()
    Y.sort()
    Z.sort()

    D = np.zeros((len(X), len(Y), len(Z)))

    for xi in range(len(X)):
        for yi in range(len(Y)):
            for zi in range(len(Z)):
                x = X[xi]
                y = Y[yi]
                z = Z[zi]
                D[xi,yi,zi] = resultsDict['[{}, {}, {}]'.format(x,y,z)]

    Dstr = np.array2string(D, max_line_width=10000000000000000, threshold=10000000000000000)

    Dstr = re.sub("\[", "{", Dstr)
    Dstr = re.sub("\]", "}", Dstr)
    Dstr = re.sub(" ", ",", Dstr)
    Dstr = re.sub("\n", ",", Dstr)

    Dstr = re.sub(",,,", ",", Dstr)
    Dstr = re.sub(",,", ",", Dstr)

    Dstr = re.sub("e-000", "", Dstr)
    Dstr = re.sub("e+000", "", Dstr)

    Dstr = re.sub("e-00", " 10^-", Dstr)
    Dstr = re.sub("e+00", " 10^", Dstr)

    Dstr = re.sub("e-0", " 10^-", Dstr)
    Dstr = re.sub("e+0", " 10^", Dstr)

    Dstr = re.sub("e-", " 10^-", Dstr)
    Dstr = re.sub("e+", " 10^", Dstr)

    with open(docname, "w") as file:
        file.write(Dstr)
    file.close()

    print("Done!")
    os.startfile(docname)
    
    
def plotFromFile(filename):
    """Plot plasma torus using saved data"""
    
    with open(filename, 'rb') as f:
        D = np.load(f)
    
    fig = mlab.figure(size=(1200, 900), bgcolor=(0,0,0))
    
    # Plot torus
    mlab.pipeline.volume( mlab.pipeline.scalar_field(D, colormap='autumn'))  # , vmin=0.3, vmax=0.6
    
    # Plot Jupiter
    plotJupiter(fig)
    mlab.show()


def main():
    """Run the code"""
    
    runProcessesModified(obj=moon, flat=flat, integrate=fieldLines, dusk=[True, False] , time=[startTime, endTime], N=N, H=H, W=W, C=C, C2min=C2min, C2max=C2max, circle=circle, phi0=phi0, theta0f=theta0f)
    
    # Or plot from file:
    #plotFromFile('plasmatorus.npy')
    
    
if __name__ == "__main__":
    main()
