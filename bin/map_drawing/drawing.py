"""
Read a cophenetic matrix to plot the distances. You can make the matrix using ete3 and tree_to_cophenetic_matrix.py
"""
import os
os.environ["CARTOPY_USER_BACKGROUNDS"] = "/home/redwards/GitHubs/crAssphage/bin/map_drawing/crassphage_maps/images"

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

def line_color(val, verbose=False):
    """
    Return the line color associated with the value.

    :param val: The numeric value, e.g. of the intensity of the line
    :param verbose: more output
    :return: the line color in hex
    """

    # current groups. This should be one less than the number of colors
    # and then we will use anything > cutoffs[-1] as the most intense color
    cutoffs = [10, 25, 50, 100]

    # current color scale. This is so we can easily change it
    colorscale = Blues

    for i, j in enumerate(cutoffs):
        if val <= j:
            return colorscale[i]
    return colorscale[-1]

def calculate_lines_dots(ll, dd, verbose=False):
    """
    Count the line data and dot data for each of our lines and dots
    :param ll: the lat lon data
    :param dd: the distances between samples
    :param verbose: more output
    :return:
    """

    # the data for our lines
    linedata = {}

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

def get_marker_size(val, verbose=False):
    """
    Get the size of the marker based on val
    :param val: the value to test
    :param verbose: more output
    :return: the size of the marker in pixels
    """

    markersizes = [4, 7, 10, 13, 16]
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
        markersize = get_marker_size(dotdata[tple], verbose)
        markercol  = line_color(dotdata[tple], verbose)

        # draw the base marker that is blck
        plt.plot(tple[0], tple[1], 'o', color='Black', markersize=markersize,
                 transform=ccrs.PlateCarree())
        if plotsingle or dotdata[tple] > 1:
            # color the line based on the value of dotat
            plt.plot(tple[0], tple[1], 'o', color=markercol, fillstyle='none', markersize=markersize+1,
                     mew=linewidth, transform=ccrs.PlateCarree())

        dotsizes[tple] = markersize

    return dotsizes

def draw_lines(linedata, plt, colorcontinents=False, plotsingle=False, verbose=False):
    """
    Draw the lines connecting the dots
    :param linedata: the line data
    :param plt: the plot image
    :param verbose: more output
    :return:
    """

    linesizes = []
    # plot the lines first so the circles are on top!
    for tple in linedata:
        if not plotsingle and linedata[tple]['count'] < 2:
            continue

        lc = line_color(linedata[tple]['count'])
        if not linedata[tple]['samecontinent'] and colorcontinents:
            sys.stderr.write("Color by continents is currently disabled. Sorry, Liz\n")

        plt.plot(linedata[tple]['x'], linedata[tple]['y'], color=lc,
                 linewidth=linedata[tple]['linewidth'],
                 zorder=linedata[tple]['count'], transform=ccrs.Geodetic())
        linesizes.append(linedata[tple]['count'])

    return linesizes


def plotmap(ll, dd, outputfile, plotsingle=False, legendfile=None, verbose=False):
    """
    Plot the map of the dna distances and lat longs
    :param ll: The lon-lats
    :param dd: The distances to use
    :param outputfile: The file name to write the image to
    :param maxdist: The maximum distance that we will scale to be maxlinewidth
    :param colorcontinents: color lines that go to different continents a different color (currently yellow)
    :param plotintensity: plot the colors of the lines by intensity vs. setting each number a color
    :param legendfile: create a separate file with the legend.
    :param linewidthbyn: scale the line width by the number of lines drawn
    :return:
    """


    if verbose:
        sys.stderr.write("Plotting the map\n")

    # there are three different alphas that we are looking at
    # the lines between samples (and potentially split those to the lines within and between continents)
    # the dots
    # the circles to themselves.

    # These are the maximum values
    #
    # Primer:     A       B      C
    # Cirlces:   530     485    289
    # Lines:      10       9     13
    # Dots:     2040    1806    680
    # at most out of these we only have 30 different numbers.

    # These numbers adjust the size of the things drawn
    ax = plt.axes(projection=ccrs.Robinson())

    # make the map global rather than have it zoom in to
    # the extents of any plotted data
    ax.set_global()
    # convert to a grayscale image. Uncomment stock_img to get color
    # ax.background_img(name='grayscale_shaded', resolution='low')
    #ax.stock_img()
    ax.coastlines()

    # figure out our dots and lines
    t = time.time()
    dotdata, linedata = calculate_lines_dots(ll, dd, verbose=verbose)
    if verbose:
        sys.stderr.write(f"Calculating the lines took {time.time() - t} seconds\n")

    t=time.time()
    dotsizes = draw_dots(dotdata, plt, plotsingle=plotsingle, verbose=verbose)
    if verbose:
        sys.stderr.write(f"Drawing the dots took {time.time() - t} seconds\n")

    linesizes  = draw_lines(linedata, plt, plotsingle=plotsingle, verbose=verbose)


    plt.savefig(outputfile)

    if legendfile:
        # create a new figure for the legend
        plt.figure(1)
        ax2 = plt.axes()
        # create the boxes for the colors

        legends = []
        labels  = []


        # combine both legends and labels to make a single legend for this figure
        alleg = legends
        allab = labels

        # ax2.legend(alleg, allab)

        markers = []
        legendtext = []
        for m in sorted(dotsizes.values()):
            emptyplot = plt.scatter([], [], s=m, marker='o', color='Black')
            markers.append(emptyplot)
            legendtext.append("")
        legendtext[0]  = "Least"
        legendtext[-1] = "Most"

        ax2.legend(markers, legendtext,
                   scatterpoints=1,
                   loc='lower left',
                   ncol=1,
                   fontsize=8)

        plt.savefig(legendfile)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot a map using ete and lat lon')
    parser.add_argument('-i', help='id.map file with lat/lon information', required=True)
    parser.add_argument('-j', help='json format of the cophenetic map file with same ids as id.map', required=True)
    parser.add_argument('-o', help='output file name', required=True)
    parser.add_argument('-a', help='alpha level for lines. Default=0.25', type=float, default=0.25)
    parser.add_argument('-l', help='linewidth for the lines connecting similar sites', default=1, type=float)
    parser.add_argument('-c', help='color the lines between continents yellow', action='store_true')
    parser.add_argument('-n', help='Plot the intensity as a fraction of max value', action='store_true')
    parser.add_argument('-g', help="Legend file. This is a simple image with the legend, and is in a separate file")
    parser.add_argument('-s', help='plot dots and lines with only a single data point (otherwise just dots)', action='store_true')
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    t = time.time()
    sys.stderr.write("Reading lat lon\n")
    lonlat = get_lon_lat(args.i, verbose=args.v)
    sys.stderr.write(f"\ttook {time.time() - t} seconds\n")

    t = time.time()
    with open(args.j, 'r') as jin:
        dist = json.load(jin)
    sys.stderr.write(f"Reading json took {time.time()-t} seconds\n")

    plotmap(lonlat, dist, args.o, plotsingle=args.s, legendfile=args.g, verbose=args.v)
