# -*- coding: utf-8 -*-
"""
Authors: Erik Hedenström and Anton Petrén, 2022
Part of the Bachelor's thesis "Using Jupiter's Moon Io as a Plasma Probe"

This file is used for plotting results.
"""

import numpy as np
import matplotlib.pyplot as plt


def singlePlot(results, display = "density"):
    """Displays the data we have collected in a 2D plot"""
    
    if display == "distance":
        dens = results[2]
        plt.xlabel("Distance from centrifugal equator (Rj)", fontsize=13)
        plt.xlim([0, 0.8])
    else:
        dens = results[0]
        plt.xlabel("Plasma density according to model ($cm^{-3}$)", fontsize=13)
        plt.xlim([1000, 2400])

    bright = results[1]
    arrayDens = np.vstack([dens, np.ones(len(dens))]).T
    arrayBright = np.array(bright)
    
    # least squares fit
    m, c = np.linalg.lstsq(arrayDens, arrayBright, rcond = None)[0]
    statistics(arrayDens, arrayBright, dens, bright, m, c)

    # make scatter plot
    #plt.scatter(dens, bright, c = results[3]) #uses colored dots to display time
    plt.scatter(dens, bright, c='C0') #no color on dots
    plt.ylim([0.3, 1.25])
    plt.plot(arrayDens, m*arrayDens + c, label='y={:.6f}x+{:.2f}'.format(m,c), c='C0')
    plt.ylabel("Observed brightness (kR)", fontsize=13)
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    plt.legend(by_label.values(), by_label.keys(), fontsize=13)
    plt.show()
    
    return


def statistics(arrayDens, arrayBright, dens, bright, m, c):
    """Print statistics"""
    
    residuals = np.linalg.lstsq(arrayDens, arrayBright, rcond = None)[1]
    lcc = np.corrcoef([dens,bright])[0][1] # linear correlation coefficient
    
    print()
    print("STATISTICS")
    print("==========")
    print("Residuals:", residuals)
    print("k =", m)
    print("c =", c)
    print("Pearson coefficient:", lcc)
    print("==========")
    
    return
