"""
Read a cophenetic matrix to plot the distances. You can make the matrix using ete3 and tree_to_cophenetic_matrix.py
"""

import os
import sys
import argparse

import gzip

import matplotlib.pyplot as plt
# set the figure size. This should be in inches?
plt.rcParams["figure.figsize"] = (22,16) # default: 22,16; for large use 88, 64
#plt.rcParams["figure.figsize"] = (88,64) # default: 22,16; for large use 88, 64
import matplotlib.colors as colors
import matplotlib.cm as cmx
from matplotlib.patches import Circle, Rectangle

import math

import cartopy.crs as ccrs
import re

country2continent = {"Albania" : "Europe", "Australia" : "Oceania", "Azerbaijan" : "Asia", "Belgium" : "Europe", "Brazil" : "South America", "Bulgaria" : "Europe", "Canada" : "North America", "China" : "Asia", "Colombia" : "South America", "Cote_dIvoire" : "Africa", "Croatia" : "Europe", "Czech_Republic" : "Europe", "Denmark" : "Europe", "Ecuador" : "South America", "Ethiopia" : "Africa", "Finland" : "Europe", "France" : "Europe", "Gambia" : "Africa", "Germany" : "Europe", "Hong_Kong" : "Asia", "Hungary" : "Europe", "Iceland" : "Europe", "India" : "Asia", "Iran" : "Asia", "Ireland" : "Europe", "Israel" : "Asia", "Italy" : "Europe", "Japan" : "Asia", "Jordan" : "Asia", "Kazakhstan" : "Asia", "Kenya" : "Africa", "Latvia" : "Europe", "Luxembourg" : "Europe", "Malaysia" : "Asia", "Malta" : "Europe", "Mexico" : "North America", "Moldova" : "Europe", "Mongolia" : "Asia", "Nepal" : "Asia", "New_Zealand" : "Oceania", "Nigeria" : "Africa", "Norway" : "Europe", "Pakistan" : "Asia", "Peru" : "South America", "Poland" : "Europe", "Portugal" : "Europe", "Russia" : "Asia", "Singapore" : "Asia", "South_Africa" : "Africa", "Spain" : "Europe", "Sri_Lanka" : "Asia", "Sudan" : "Africa", "Sweden" : "Europe", "Switzerland" : "Europe", "The_Netherlands" : "Europe", "United_States" : "North America", "USA" : "North America", "Zambia" : "Africa", "crAssphage_A" : "crAssphage_A"}


def get_lon_lat(idf, maxtoget=50000):
    """
    Get the longitude and latitude of different ids. Note that we have longitude first to work with cartopy
    :param idf: the id.map file
    :param maxtoget: the maxiumum number of ids to get. This is just for debugging
    :return:
    """
    lonlat = {}
    count = 0
    global verbose
    with open(idf, 'r') as fin:
        for l in fin:
            if count > maxtoget:
                break
            count+=1
            s=re.search('latitude=(\S+)\]', l)
            if not s:
                sys.stderr.write("No latitude in {}".format(l))
                continue
            lat=s.group(1)

            s = re.search('longitude=(\S+)\]', l)
            if not s:
                sys.stderr.write("No longitude in {}".format(l))
                continue
            lon = s.group(1)
            p=l.split("\t")

            try:
                lat = float(lat)
                lon = float(lon)
            except:
                sys.stderr.write("There was an error parsing the latitude and longitude from {}\n".format(l))
                continue

            lonlat[p[0]] = (lon, lat)
            newname = p[0].replace('|', '_')
            lonlat[newname] = (lon, lat)
    return lonlat

def latlon2distance(lat1, long1, lat2, long2, miles=False):
    """Convert two coordinates to distance.

    This is an approximation since the earth is not spherical, but accuracy is <100m, especially for close points

    This code was taken from http://www.johndcook.com/python_longitude_latitude.html

    Latitude is measured in degrees north of the equator; southern locations have negative latitude.
    Similarly, longitude is measured in degrees east of the Prime Meridian. A location 10deg west of
    the Prime Meridian, for example, could be expressed as either 350deg  east or as -10deg east.

    Arguments: lat1, long1; lat2, long2; miles is a boolean. If you want miles set it to true. Else set it to false

    """
    global verbose

    if lat1 == lat2 and long1 == long2:
        return 0


    # Convert latitude and longitude to
    # spherical coordinates in radians.
    degrees_to_radians = math.pi / 180.0

    # phi = 90 - latitude
    phi1 = (90.0 - lat1) * degrees_to_radians
    phi2 = (90.0 - lat2) * degrees_to_radians

    # theta = longitude
    theta1 = long1 * degrees_to_radians
    theta2 = long2 * degrees_to_radians

    # Compute spherical distance from spherical coordinates.

    # For two locations in spherical coordinates
    # (1, theta, phi) and (1, theta, phi)
    # cosine( arc length ) =
    #    sin phi sin phi' cos(theta-theta') + cos phi cos phi'
    # distance = rho * arc length

    cos = (math.sin(phi1) * math.sin(phi2) * math.cos(theta1 - theta2) + math.cos(phi1) * math.cos(phi2))
    try:
        arc = math.acos(cos)
    except Exception as err:
        sys.stderr.write("There was an err: {} trying to take the acos of ({})\n".format(err, cos))
        arc=0
    # Remember to multiply arc by the radius of the earth
    # in your favorite set of units to get length.
    #
    # To convert to miles multiple arc by 3960
    # To convert to kilometers multiply arc by 6373

    if miles:
        arc *= 3960
    else:
        arc *= 6373

    return arc


def closest_dna_dist(matrixfile):
    """
    Read the matrix file and get the id of the point with the closest distance that is not ourself
    :param treefile: The cophenetic matrix file to read
    :return: a dict of a node and its closest leaf
    """

    global verbose
    if verbose:
        sys.stderr.write("Getting closest distances\n")
    distances = {}

    if matrixfile.endswith('.gz'):
        with gzip.open(matrixfile, 'rt') as f:
            l = f.readline()
            ids = l.rstrip().split("\t")
            for i,name in enumerate(ids):
                if i == 0:
                    continue
                distances[name] = {}
            for l in f:
                data = l.rstrip().split("\t")
                for i,dist in enumerate(data):
                    if i == 0:
                        continue
                    distances[data[0]][ids[i]] = float(dist)
                    distances[ids[i]][data[0]] = float(dist)
    else:
        with open(matrixfile, 'r') as f:
            l = f.readline()
            ids = l.rstrip().split("\t")
            for i,name in enumerate(ids):
                if i == 0:
                    continue
                distances[name] = {}
            for l in f:
                data = l.rstrip().split("\t")
                for i,dist in enumerate(data):
                    if i == 0:
                        continue
                    distances[data[0]][ids[i]] = float(dist)
                    distances[ids[i]][data[0]] = float(dist)


    closest = {}
    for d in distances:
        closest[d] = {}
        for k in sorted(distances[d], key=distances[d].get):
            if k == d:
                continue
            closest[d][k] = distances[d][k]
            break
    if verbose:
        sys.stderr.write("From\tTo\tDistance\n")
        for d in distances:
            for k in closest[d]:
                sys.stderr.write("{}\t{}\t{}\n".format(d, k, closest[d][k]))


    if verbose:
        sys.stderr.write("\n\n\nDone\n")
    return closest

def plotmap(ll, dd, outputfile, alpha, linewidth=1, bounds=None, maxdist=1, maxlinewidth=6, colorcontinents=False):
    """
    Plot the map of the dna distances and lat longs
    :param ll: The lon-lats
    :param dd: The distances to use
    :param outputfile: The file name to write the image to
    :param bounds: the boundary for the lat long as a 4-ple array
    :param maxdist: The maximum distance that we will scale to be maxlinewidth
    :param colorcontinents: color lines that go to different continents a different color (currently yellow)
    :return:
    """
    global verbose

    if verbose:
        sys.stderr.write("Plotting the map\n")

    ax = plt.axes(projection=ccrs.Robinson())

    # make the map global rather than have it zoom in to
    # the extents of any plotted data
    if not bounds:
        ax.set_global()

    ax.stock_img()
    ax.coastlines()

    ## color the lines based on the maximum distance value
    jet = cm = plt.get_cmap('jet')
    cNorm = colors.Normalize(vmin=0, vmax=maxdist)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)

    # Using contourf to provide my colorbar info, then clearing the figure
    Z = [[0, 0], [0, 0]]
    levels = range(0, int(100 * maxdist) + 10, 10)
    CS3 = plt.contourf(Z, levels, cmap=jet)
#    plt.clf()


    # NOTE: longitude before latitude!!
    # plt.plot([sdlon, brislon], [sdlat, brislat], color='blue', linewidth=2,  transform=ccrs.Geodetic())

    # plot the circles for each sample site
    # markerfacecolor="None",

    for lonlat in ll.values():
        if bounds and ((lonlat[1] < bounds[0] or lonlat[1] > bounds[2]) or (lonlat[0] < bounds[1] or lonlat[0] > bounds[3])):
            if verbose:
                sys.stderr.write("Not in bounding box: {}\n".format(lonlat))
            continue
        if verbose:
            sys.stderr.write("Kept location: {}\n".format(lonlat))
        plt.plot(lonlat[0], lonlat[1], 'o', color='Black', alpha=alpha, markersize=10, transform=ccrs.PlateCarree())

    # how many lines and circles do we draw?
    circleat = {}
    lineat   = {}
    for idx1 in dd:
        for idx2 in dd[idx1]:
            # this should only happen when we do best DNA distances
            if idx1 not in ll:
                sys.stderr.write("NO Lat/Lon for {}\n".format(idx1))
                continue
            if idx2 not in ll:
                sys.stderr.write("NO Lat/Lon for {}\n".format(idx2))
                continue

            linecolor = 'red'
            scaledalpha = alpha
            if colorcontinents:
                # figure out if they are from the same continent
                m = re.search('\d{8}_(\w+)\_\d', idx1)
                cont1 = country2continent.get(m.groups(0)[0], "unknown")
                m = re.search('\d{8}_(\w+)\_\d', idx2)
                cont2 = country2continent.get(m.groups(0)[0], "unknown")
                if cont1 != cont2:
                    linecolor = 'yellow'
                    scaledalpha = alpha * 0.75

            if bounds and ((ll[idx1][1] < bounds[0] or ll[idx1][1] > bounds[2]) or (ll[idx1][0] < bounds[1] or ll[idx1][0] > bounds[3])):
                if verbose:
                    sys.stderr.write("{} out of bounds. Skipped\n".format(idx1))
                continue

            if bounds and ((ll[idx2][1] < bounds[0] or ll[idx2][1] > bounds[2]) or (ll[idx2][0] < bounds[1] or ll[idx2][0] > bounds[3])):
                if verbose:
                    sys.stderr.write("{} out of bounds. Skipped\n".format(idx2))
                continue

            if linewidth == 0:
                linewidth = dd[idx1][idx2]
                linewidth = (linewidth/maxdist) * maxlinewidth

            if verbose:
                sys.stderr.write("{} to {}: distance: {} km. Genetic distance {}. Line width {}\n".format(
                    idx1, idx2, latlon2distance(ll[idx1][1], ll[idx1][0], ll[idx2][1], ll[idx2][0]), dd[idx1][idx2], linewidth))


            #colorVal = scalarMap.to_rgba(dd[idx1][idx2])

            if latlon2distance(ll[idx1][1], ll[idx1][0], ll[idx2][1], ll[idx2][0]) < 100:
                if verbose:
                    sys.stderr.write("Adding a circle for {} and {}\n".format(ll[idx1][0], ll[idx1][1]))
                # add a red circle for this object.
                # we need to use some simple trig to find the center of the circle whose point on the circumference
                # is at our lat lon
                radius = 3
                if bounds:
                    radius = 1.5
                circlon = ll[idx1][0] - (radius * math.sin(2 * math.pi))
                circlat = ll[idx1][1] - (radius * math.cos(2 * math.pi))

                circ = Circle((circlon, circlat), transform=ccrs.Geodetic(), radius=radius,
                                     linewidth=linewidth, alpha=scaledalpha, color=linecolor, fill=False)
                circleat[(circlon, circlat)] = circleat.get((circlon, circlat), 0) + 1
                ax.add_artist(circ)
            else:
                # plot a red line between two points
                plt.plot([ll[idx1][0], ll[idx2][0]], [ll[idx1][1], ll[idx2][1]], color=linecolor, linewidth=linewidth,
                         alpha=scaledalpha, transform=ccrs.Geodetic())
                lineat[((ll[idx1][0], ll[idx2][0]), (ll[idx1][1], ll[idx2][1]))] = lineat.get(((ll[idx1][0], ll[idx2][0]), (ll[idx1][1], ll[idx2][1])), 0) + 1

    #    plt.colorbar(CS3)

    # add a color bar for the lines and circles
    rect = Rectangle((10, 10), 140, 130, linewidth=5, edgecolor='b', facecolor='none')
    #plt.legend(handles=[rect])

    #grad = plt.imshow([[0.,1.], [0.,1.]], cmap = plt.cm.Reds, interpolation = 'bicubic')
    #plt.legend(handles=[grad])
    fig = plt.figure()
    ax2 = fig.add_axes()
    ax2.add_patch(rect)


    #plt.show()
    plt.savefig(outputfile)

    # sys.stderr.write("We drew a max of {} circles\n".format(max(circleat.values())))
    # sys.stderr.write("And we drew a max of {} lines\n".format(max(lineat.values())))
    sys.stderr.write("Circles,{}\nLines,{}".format(",".join(map(str, circleat.values())), ",".join(map(str, lineat.values()))))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot a map using ete and lat lon')
    parser.add_argument('-i', help='id.map file with lat/lon information', required=True)
    parser.add_argument('-m', help='cophenetic map file with same ids as id.map', required=True)
    parser.add_argument('-o', help='output file name', required=True)
    parser.add_argument('-a', help='alpha level for lines. Default=0.25', type=float, default=0.25)
    parser.add_argument('-l', help='linewidth for the lines connecting similar sites', default=1, type=float)
    parser.add_argument('-c', help='color the lines between continents yellow', action='store_true')
    parser.add_argument('-b', help='geographic bounds. Use top left, bottom right. e.g. 75,35:35,-25')
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    global verbose
    verbose = False
    if args.v:
        verbose = True

    bounds = None
    if args.b:
        ll = args.b.split(":")
        topleft = ll[0].split(",")
        bottomright = ll[1].split(",")
        try:
            bounds=[int(topleft[0]), int(topleft[1]), int(bottomright[0]), int(bottomright[1])]
            if bounds[0] > bounds[2]:
                (bounds[0], bounds[2])=(bounds[2], bounds[0])
            if bounds[1] > bounds[3]:
                (bounds[1], bounds[3])=(bounds[3], bounds[1])
        except:
            sys.stderr.write("There was an error parsing integers from {}. Please do not include N/E/S/W, just +/- ints\n".format(args.b))
            sys.exit()

    lonlat = get_lon_lat(args.i)
    # dist = best_dna_dist(get_dna_distance(args.t))
    dist = closest_dna_dist(args.m)
    plotmap(lonlat, dist, args.o, args.a, linewidth=args.l, bounds=bounds, colorcontinents=args.c)