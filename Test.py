# -*- coding: utf-8 -*-
"""
Authors: Erik Hedenström and Anton Petrén, 2022
Part of the Bachelor's thesis "Using Jupiter's Moon Io as a Plasma Probe"

This file is used to run the code. Choose desired parameters and import 
observations of aurora on Io.
"""

import ObjectPosition
import PlasmaTorus
import PlotResults
import multiprocessing as mp
from datetime import datetime

# Import txt files with aurora observations on Io. Each line in the txt file
# should be in the format
# "o49d01010   1997-10-14T02:52:32.52       353.60992       241.30518     0.617715"
filenamesIo = ["data/KEXvisit3_s.txt", "data/KEXvisit4_s.txt", "data/KEXvisit6_s.txt", "data/KEXvisit7_s.txt", "data/KEXvisit8_s.txt", "data/KEXvisit9_s.txt", "data/KEXvisit10_s.txt", "data/KEXvisit12_s.txt", "data/KEXvisit13_s.txt", "data/KEXvisit14_s.txt", "data/KEXvisit15_s.txt", "data/KEXvisit16_s.txt", "data/KEXvisit17_s.txt", "data/KEXvisit22_s.txt", "data/KEXvisit24_s.txt", "data/KEXvisit25_s.txt", "data/KEXvisit26_s.txt", "data/KEXvisit27_s.txt"]

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
C2min = 5.6  # Ribbon minimum distance from center of Jupiter. Overrides C[1]
C2max = 5.6  # Ribbon maximum distance from center of Jupiter. Overrides C[1]
phi0 = 0  # Deviation from true dawn/dusk-direction (degrees)
# It is assumed that C2max is on the dusk side and that C2min is on the dawn
# side. The value provided in C[1] above should be the average ribbon distance,
# smaller than C2max and larger than C2min. All other C values will be updated
# such that the whole torus moves with the position of the ribbon.
# If C2min=C2max, dawn-dusk effects will be ignored and the C[1] value above
# will be used for the ribbon.

# Data selection
dusk = [True, False] 
# [True] = plot only dusk data
# [False] = plot only dawn data
# [True, False] = plot dusk and dawn data

# Plotting options
x_axis = ""  # Empty string = plasma density, "distance" = distance from centrifugal equator
plotresult = True  # Plot the results

# SPICE options
startTime = datetime(year=1980, month=8, day=23, hour=19, minute=15, second=0)
endTime = datetime(year=2022, month=8, day=24, hour=1, minute=0, second=0)
moon = "Io"


def runProcesses(obj, flat, integrate, dusk, time, N, H, W, C, C2min, C2max, circle, phi0, theta0f):
    """This main function serves only as a function to start processes and graph the data they return, all calculation are done in other pocesses"""
    
    with mp.Manager() as manager:
        processes = []
        lock = mp.Lock()

        if obj == "Io":
            data = unpackDataIo()
        else:
            data = unpackDataEuropa()

        results = manager.list([])

        for i in range(NUMBEROFPROCESSES):
            p = mp.Process(target = testChain, args = (obj, data, lock, i, NUMBEROFPROCESSES, results, flat, integrate, N, H, W, C, C2min, C2max, circle, phi0, theta0f))
            processes.append(p)

        for p in processes:
            p.start()

        for p in processes:
            p.join()

        resultsDict = dict(results)

        brightness = []
        density = []
        height = []
        resultsData = []
        seconds = []
        for i in range(len(results)):
            if resultsDict[i][4] in dusk:
                timedelta1 = (resultsDict[i][1]-time[0]).total_seconds()
                timedelta2 = (resultsDict[i][1]-time[1]).total_seconds()
                if timedelta1>0 and timedelta2<0:
                    resultsData.append(resultsDict[i])
                    brightness.append(resultsDict[i][2])
                    density.append(resultsDict[i][3])
                    height.append(resultsDict[i][5])
                    seconds.append((resultsDict[i][1]-time[0]).total_seconds())

        return [density, brightness, height, seconds, resultsData]


def testChain(obj, data, lock, i, n, results, flat, integrate, N, H, W, C, C2min, C2max, circle, phi0, theta0f):
    """Allows us to use a single process to run multiple tests, lowering the minimum number of processes"""
    
    for x in range(i, len(data), n):
        result = singleTest(obj, data[x], flat, integrate, N, H, W, C, C2min, C2max, circle, phi0, theta0f)
        lock.acquire()
        results.append((x, result))
        lock.release()


def singleTest(obj, dataPoint, flat, integrate, N, H, W, C, C2min, C2max, circle, phi0, theta0f):
    """Runs the code we wish to use to test a single data point and returns the value we wish to display"""

    dateStr = dataPoint[1]
    y = int(dateStr[:4])
    m = int(dateStr[5:7])
    day = int(dateStr[8:10])
    h = int(dateStr[11:13])
    min = int(dateStr[14:16])
    s = int(float(dateStr[17:]))
    d = datetime(year=y, month=m, day=day, hour=h, minute=min, second=s)
    dataPoint = [dataPoint[0], d, float(dataPoint[2])*(-1), float(dataPoint[3])*(-1), float(dataPoint[4])]

    name = dataPoint[0]  # name of observation
    d = ObjectPosition.localTime("Jupiter", dataPoint[1])  # corresponding jupiter time
    pos = ObjectPosition.Spice(obj, d)
    posPolar = ObjectPosition.cartesianToPolar(pos)
    sunpos = ObjectPosition.Spice("Sun", d)
    phiSun = ObjectPosition.cartesianToPolar(sunpos)[2]
    arg = [posPolar[0], posPolar[1], posPolar[2], phiSun]
    density, eqDev = PlasmaTorus.plasmaDensity(coordinates=arg, flat=flat, integrate=integrate, N=N, H=H, W=W, C=C, C2min=C2min, C2max=C2max, circle=circle, phi0=phi0, theta0f=theta0f)
    phi1 = posPolar[2]
    phi2 = 0  # TODO: compute phi2 (for verification)
    spiceEv = [name, d, phi1, phi2, density]

    op = [name, d, dataPoint[4], density, PlasmaTorus.isDusk(phi1, phiSun, phi0), abs(eqDev)]

    return op


def unpackDataIo():
    """Reads data from the data files and places it in an array"""
    
    stringlist = []
    for file in filenamesIo:
        opened = open(file)
        openedLines = opened.readlines()
        stringlist += openedLines
        opened.close()

    data = []
    for string in stringlist:
        data.append(stringToList(string))

    return data


def stringToList(string):
    """Convert a string to a list"""
    
    tag = string[0:9]
    date = string[12:34]
    theta = string[41:51]
    earthTheta = string[57:66]
    brightness = string[71:79]
    return [tag, date, theta, earthTheta, brightness]


def main():
    """Run the code"""
    
    data1 = runProcesses(obj=moon, flat=flat, integrate=fieldLines, dusk=dusk, time=[startTime, endTime], N=N, H=H, W=W, C=C, C2min=C2min, C2max=C2max, circle=circle, phi0=phi0, theta0f=theta0f)
    if plotresult:
        PlotResults.singlePlot(data1, x_axis)
        #print(data1)
    
    
if __name__ == "__main__":
    main()
