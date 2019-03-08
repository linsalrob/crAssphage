"""
Read a cophenetic matrix to plot the distances as pie charts
"""
import os
os.environ["CARTOPY_USER_BACKGROUNDS"] = "/home/redwards/GitHubs/crAssphage/bin/map_drawing/crassphage_maps/images"

import sys
import argparse
import matplotlib.pyplot as plt
# set the figure size. This should be in inches?
#plt.rcParams["figure.figsize"] = (22,16) # default: 22,16; for large use 88, 64

#plt.rcParams["figure.figsize"] = (44,32) # default: 22,16; for large use 88, 64

plt.rcParams["figure.figsize"] = (88,64) # default: 22,16; for large use 88, 64
import matplotlib.colors as colors
import matplotlib.cm as cmx
from matplotlib.patches import Circle, Rectangle
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

import math
import json
import time

import pprint

import cartopy.crs as ccrs
import re

import numpy as np
from itertools import compress

from crassphage_maps import closest_dna_dist, get_lon_lat, latlon2distance, country2continent
from crassphage_maps import green2red, green2yellow, evenly_select, GnBu5, Blues, YlOrBr, YlOrRd
from roblib import bcolors

def cols(val, verbose=False, getscale=False):
    """
    Return the line color associated with the value.

    :param val: The numeric value, e.g. of the intensity of the line
    :param verbose: more output
    :param getscale: just return the color scale array, regardless of val
    :return: the line color in hex
    """

    # current groups. This should be one less than the number of colors
    # and then we will use anything > cutoffs[-1] as the most intense color
    cutoffs = [10, 25, 50, 100]

    # current color scale. This is so we can easily change it
    colorscale = YlOrRd

    if getscale:
        return colorscale

    for i, j in enumerate(cutoffs):
        if val <= j:
            return colorscale[i]
    return colorscale[-1]


def get_pie_size(val, verbose=False, getscale=False):
    """
    Get the size of the pie chart based on val.
    :param val: the value to test
    :param verbose: more output
    :param getscale: return the marker sizes regardless of val
    :return: the size of the marker in pixels
    """

    # these are the sizes of the pies
    # piesizes = [22, 25, 28, 31, 34]
    # piesizes = [15, 18, 21, 24, 27]
    # piesizes = [4, 9, 14, 19, 25]
    # piesizes = [13, 16, 19, 22, 25]
    # piesizes = [11, 14, 17, 20, 23]
    # piesizes = [9, 12, 15, 18, 21]
    # piesizes = [8, 11, 14, 17, 20]
    # piesizes = [6, 9, 12, 15, 18]
    # piesizes = [4, 6, 8, 10, 14]
    # piesizes = [2, 4, 6, 8, 12]
    # piesizes = [0.5, 1, 1.5, 2, 2.5]
    piesizes = [0.5, 0.7, .9, 1.1, 1.3]
    # there should be one less maxpievals than piesizes and then we use piesizes[-1] for anything larger
    maxpievals = [10, 20, 30, 40]

    if getscale:
        return piesizes

    for i,m in enumerate(maxpievals):
        if val <= m:
            return piesizes[i]
    return piesizes[-1]


def get_alpha(val, verbose=False, getscale=False):
    """
    Get an alpha level associated with this amount
    :param val: the value
    :param verbose: more output
    :param getscale: get the alpha scale regardless of val
    :return: the alpha level, a number between 0 and 1
    """

    alphalevels = [0.2, 0.4, 0.6, 0.8, 1]
    if verbose:
        sys.stderr.write(f"Getting alpha total: {val}\n")
    alphavals   = [10,20,30,40]

    if getscale:
        return alphalevels

    for i, m in enumerate(alphavals):
        if val <= m:
            return alphalevels[i]
    return alphalevels[-1]


def calculate_samediff(ll, dd, samekm=100, verbose=False):
    """
    Calculate same/different values for each location. We also concatenate all the locations
    < samekm apart.
    :param ll: the lat lon data
    :param dd: the distances between samples
    :param samekm: the distance between two sites to call them the same
    :param verbose: more output
    :return:
    """

    data = {}
    # first we want to make sure we have latitude and longitude for everything
    missingll = False
    for idx1 in dd:
        if idx1 not in ll:
            if verbose:
                sys.stderr.write("NO Lat/Lon for {}\n".format(idx1))
            missingll = True
    if missingll:
        sys.stderr.write(f"{bcolors.RED}FATAL: We are missing latitudes and longitudes. You should check this{bcolors.ENDC}\n")
        sys.exit(2)

    # first we concatenate everything
    # our data structure is a hash of locations : {location : distance}
    # each location has a single value that is the single closest location and its genetic distance
    # first we check if they are the same or different locations, then we check all locations that we've
    # already counted to see if they are within samekm km. If so we add the same/diff to that. Otherwise we
    # set it as a new location

    for idx1 in dd:
        ll1 = (ll[idx1][0], ll[idx1][1])
        # we set two values as 1/0 to increment either in the dict as appropriate
        sameloc = 0
        diffloc = 0
        idx2 = list(dd[idx1].keys())[0] # there is only one key in the dict!
        ll2 = (ll[idx2][0], ll[idx2][1])
        if latlon2distance(ll[idx1][1], ll[idx1][0], ll[idx2][1], ll[idx2][0]) < samekm:
            sameloc = 1
            diffloc = 0
        else:
            sameloc = 0
            diffloc = 1
        added = False
        for s in data:
            if latlon2distance(ll[idx1][1], ll[idx1][0], data[s]['lon'], data[s]['lat']) < samekm:
                added = True
                data[s]['same'] += sameloc
                data[s]['diff'] += diffloc
                break
        if not added:
            data[idx1] = {
                'lon': ll[idx1][1],
                'lat': ll[idx1][0],
                'same': sameloc,
                'diff': diffloc
            }

    # figure out the maximum total
    # and add totals to everything
    maxt=0
    for s in data:
        data[s]['total'] = data[s]['same'] + data[s]['diff']
        if data[s]['total'] > maxt:
            maxt = data[s]['total']

    if verbose:
        p = pprint
        sys.stderr.write(f"SAMEDIFF\n{p.pformat(data)}\n")

    return data


def plot_pie_inset(data,lat,lon,ax,width,c):
    """
    Plot the pie charts on the axis.
    This is taken from
    https://stackoverflow.com/questions/45266955/adding-pie-chart-at-given-coordinates-to-cartopy-projection
    :param data: The array of the numbers for the pie chart
    :param lon: longitude
    :param lat: latitude
    :param ax: axes to plot upon
    :param width: size of the pie. we have to set height/width to this
    :param c: colors for the wedges
    :return:
    """
    lonr, latr = ccrs.Robinson().transform_point(lon, lat, ccrs.PlateCarree())
    ax_sub= inset_axes(ax, width=width, height=width, loc=10,
                       bbox_to_anchor=(lonr, latr),
                       bbox_transform=ax.transData,
                       borderpad=0)
    wedges,texts= ax_sub.pie(data, startangle=90, colors=c)
    # set some transparancy
    #al = get_alpha(sum(data))
    al = 0.5
    for i in range(len(wedges)):
        wedges[i].set_alpha(al)

    ax_sub.set_aspect("equal")

def add_pies(data, x, verbose=False):
    """
    Add the pie charts
    :param data: our same/difference/total/lat/lon data
    :param x: the plot axes/fram
    :param verbose: more output
    :return:
    """

    for d in data:
        o = data[d]['lon']
        a = data[d]['lat']
        t = data[d]['total']
        s = data[d]['same']
        f = data[d]['diff']

        c = cols(0, getscale=True)

        colpos = (0,-1)
        if 0 == s and 1 ==f:
            colpos = (len(c) // 2, len(c) // 2)


        z = get_pie_size(t)
        if verbose:
            sys.stderr.write(f"Total: {t} Size of pie: {z}\n")
        plot_pie_inset([f, s], o, a, x, z, [c[colpos[0]], c[colpos[-1]]])


def plotmap(ll, dd, outputfile, plotsingle=False, dotalpha=False, legendfile=None, verbose=False):
    """
    Plot the map of the dna distances and lat longs
    :param ll: The lon-lats
    :param dd: The distances to use
    :param outputfile: The file name to write the image to
    :param maxdist: The maximum distance that we will scale to be maxlinewidth
    :param colorcontinents: color lines that go to different continents a different color (currently yellow)
    :param dotalpha: use an alpha on the dots
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
    #ax.background_img(name='grayscale_shaded', resolution='low')
    ax.background_img(name='grey_blue_ocean6', resolution='low')

    #ax.stock_img()
    #ax.coastlines()

    # figure out our dots and lines
    t = time.time()
    data = calculate_samediff(ll, dd, samekm=250, verbose=verbose)
    if verbose:
        sys.stderr.write(f"Calculating the lots took {time.time() - t} seconds\n")

    t=time.time()
    add_pies(data, ax, verbose=verbose)

    makelegend = True
    if makelegend:
        ll = -10
        ldata = {'one' : {'lon': -2, 'lat': 0, 'same': 0, 'diff': 1, 'total': 1}}
        for ps in [2, 6, 12, 16, 24]:
            ldata[ps] = {
                'lon': ll,
                'lat': 0,
                'same': ps,
                'diff': ps,
                'total': 2*ps
            }
            ll -= 12
        add_pies(ldata, ax, verbose=verbose)

    if verbose:
        sys.stderr.write(f"Drawing the pies took {time.time() - t} seconds\n")

    plt.savefig(outputfile)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot a map using ete and lat lon')
    parser.add_argument('-i', help='id.map file with lat/lon information', required=True)
    parser.add_argument('-j', help='json format of the cophenetic map file with same ids as id.map', required=True)
    parser.add_argument('-o', help='output file name', required=True)
    parser.add_argument('-a', help='use an alpha for dots. Default=False', action='store_true')
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

    plotmap(lonlat, dist, args.o, plotsingle=args.s, legendfile=args.g, dotalpha=args.a, verbose=args.v)
