"""
In this example, we only have dots. The Size of the dots are the number of strains at each point
and the color of the dots are the number of connnections to that dot.
"""

import os
import sys
import argparse

import matplotlib.pyplot as plt
# set the figure size. This should be in inches?
plt.rcParams["figure.figsize"] = (22,16) # default: 22,16; for large use 88, 64
#plt.rcParams["figure.figsize"] = (88,64) # default: 22,16; for large use 88, 64
import matplotlib.colors as colors
import matplotlib.cm as cmx
from matplotlib.patches import Circle, Rectangle

import math
import json
import time

import cartopy.crs as ccrs
import re

import numpy as np
from itertools import compress

from crassphage_maps import closest_dna_dist, get_lon_lat, latlon2distance, country2continent
from crassphage_maps import green2red, green2yellow, evenly_select, GnBu5, Blues
from roblib import bcolors


def calculate_lines_dots(ll, dd, verbose=False):
    """
    Count the line data and dot data for each of our lines and dots
    :param ll: the lat lon data
    :param dd: the distances between samples
    :param verbose: more output
    :return: dotdata: dict of counts of that dot. Dotsim = counts of incoming connections to that dot (including
    from itself)
    """

    # the data for our lines
    dotsim = {}

    # first we want to make sure we have latitude and longitude for everything
    missingll = False
    for idx1 in dd:
        if idx1 not in ll:
            if verbose:
                sys.stderr.write("NO Lat/Lon for {}\n".format(idx1))
            missingll = True
    if missingll:
        sys.stderr.write(f"{bcolors.FATAL}FATAL: We are missing latitudes and longitudes. You should check this{bcolors.ENDC}\n")
        sys.exit(2)

    # note that now we calculate where everything should be and then plot it based on maximum values!
    dotdata = {}
    seen = set()
    for idx1 in dd:
        ll1 = (ll[idx1][0], ll[idx1][1])
        if ll1 not in dotdata:
            dotdata[ll1] = 0
        for idx2 in dd[idx1]:
            dotdata[ll1] += 1
            ll2 = (ll[idx2][0], ll[idx2][1])
            # check we only see each pair once
            linecoordsAB = "\t".join(map(str, [ll[idx1][0], ll[idx2][0], ll[idx1][1], ll[idx2][1]]))
            linecoordsBA = "\t".join(map(str, [ll[idx1][1], ll[idx2][1], ll[idx1][0], ll[idx2][0]]))
            #if linecoordsAB in seen or linecoordsBA in seen:
            #    continue
            seen.add(linecoordsAB)
            seen.add(linecoordsBA)

            # if the distance between them is less than 100km merge to a single point, and so we
            # just use idx1
            if latlon2distance(ll[idx1][1], ll[idx1][0], ll[idx2][1], ll[idx2][0]) < 100:
                # are considered the same location, so we increment dotat[idx1]
                dotdata[ll1] += 1
                continue

            samecontinent = True
            # figure out if they are from the same continent
            m = re.search(r'\d{8}_(\w+)_\d', idx1)
            cont1 = country2continent.get(m.groups(0)[0], "unknown")
            m = re.search(r'\d{8}_(\w+)_\d', idx2)
            cont2 = country2continent.get(m.groups(0)[0], "unknown")
            if cont1 != cont2:
                samecontinent = False

            if linecoordsAB in linedata:
                linedata[linecoordsAB]['count'] +=1
            else:
                linedata[linecoordsAB] = {
                    'count' : 1,
                    'x': [ll[idx1][0], ll[idx2][0]],
                    'y': [ll[idx1][1], ll[idx2][1]],
                    'samecontinent': samecontinent,
                    'linewidth' : 2
                }
    if verbose:
        sys.stderr.write("Dots,{}\n".format(",".join(map(str, dotdata.values()))))
        l = [linedata[x]['count'] for x in linedata]
        sys.stderr.write("Lines:{}\n".format(",".join(map(str, l))))
    return dotdata, linedata


def get_dot_size(val, verbose=False):
    """
    Get the size of the marker based on val
    :param val: the value to test
    :param verbose: more output
    :return: the size of the marker in pixels
    """

    markersizes = [6, 9, 12, 15, 18]
    # there should be one less maxmakervals than markersizes and then we use markersizes[-1] for anything larger
    maxmarkervals = [10, 20, 30, 40]

    for i,m in enumerate(maxmarkervals):
        if val <= m:
            return markersizes[i]
    return markersizes[-1]



def draw_dots(dotdata, plt, linewidth=2, plotsingle=False, verbose=False):
    """
    draw just the dots
    :param dotdata: the counts of occurrences of the dots
    :param plt: the matplotlib plot
    :param plotsingle: plot dots with only a single sequence
    :param verbose: more output
    :return: the datalegend and the datalabels
    """

    dotlegend = []
    dotlabels = []
    dotsizes = {}


    for tple in sorted(dotdata, key=dotdata.get):
        markersize = get_dot_size(dotdata[tple], verbose)
        markercol  = get_dot_color(dotdata[tple], verbose)

        # draw the base marker that is blck
        plt.plot(tple[0], tple[1], 'o', color='Black', markersize=markersize,
                 transform=ccrs.PlateCarree())
        if plotsingle or dotdata[tple] > 1:
            # color the line based on the value of dotat
            plt.plot(tple[0], tple[1], 'o', color=markercol, fillstyle='none', markersize=markersize+1,
                     mew=linewidth, transform=ccrs.PlateCarree())

        dotsizes[tple] = markersize

    return dotsizes




if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-f', help='', required=True)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()
