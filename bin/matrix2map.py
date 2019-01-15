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

import numpy as np
from itertools import compress

country2continent = {"Albania": "Europe", "Australia": "Oceania", "Azerbaijan": "Asia", "Belgium": "Europe",
                     "Brazil": "South America", "Bulgaria": "Europe", "Canada": "North America", "China": "Asia",
                     "Colombia": "South America", "Cote_dIvoire": "Africa", "Croatia": "Europe",
                     "Czech_Republic": "Europe", "Denmark": "Europe", "Ecuador": "South America",
                     "Ethiopia": "Africa", "Finland": "Europe", "France": "Europe", "Gambia": "Africa",
                     "Germany": "Europe", "Hong_Kong": "Asia", "Hungary": "Europe", "Iceland": "Europe",
                     "India": "Asia", "Iran": "Asia", "Ireland": "Europe", "Israel": "Asia", "Italy": "Europe",
                     "Japan": "Asia", "Jordan": "Asia", "Kazakhstan": "Asia", "Kenya": "Africa", "Latvia": "Europe",
                     "Luxembourg": "Europe", "Malaysia": "Asia", "Malta": "Europe", "Mexico": "North America",
                     "Moldova": "Europe", "Mongolia": "Asia", "Nepal": "Asia", "New_Zealand": "Oceania",
                     "Nigeria": "Africa", "Norway": "Europe", "Pakistan": "Asia", "Peru": "South America",
                     "Poland": "Europe", "Portugal": "Europe", "Russia": "Asia", "Singapore": "Asia",
                     "South_Africa": "Africa", "Spain": "Europe", "Sri_Lanka": "Asia", "Sudan": "Africa",
                     "Sweden": "Europe", "Switzerland": "Europe", "The_Netherlands": "Europe",
                     "United_States": "North America", "USA": "North America", "Zambia": "Africa",
                     "crAssphage_A": "crAssphage_A"}

# these are arrays of 50 colors from http://www.perbang.dk/rgbgradient/ between FF0000 (red) and 0000FF (blue) or
# #00FF00 (green) and #FFFF00 (yellow)
green2yellow = ["#00FF00", "#05FF00", "#0AFF00", "#0FFF00", "#14FF00", "#1AFF00", "#1FFF00", "#24FF00", "#29FF00",
                "#2EFF00", "#34FF00", "#39FF00", "#3EFF00", "#43FF00", "#48FF00", "#4EFF00", "#53FF00", "#58FF00",
                "#5DFF00", "#62FF00", "#68FF00", "#6DFF00", "#72FF00", "#77FF00", "#7CFF00", "#82FF00", "#87FF00",
                "#8CFF00", "#91FF00", "#96FF00", "#9CFF00", "#A1FF00", "#A6FF00", "#ABFF00", "#B0FF00", "#B6FF00",
                "#BBFF00", "#C0FF00", "#C5FF00", "#CAFF00", "#D0FF00", "#D5FF00", "#DAFF00", "#DFFF00", "#E4FF00",
                "#EAFF00", "#EFFF00", "#F4FF00", "#F9FF00", "#FFFF00"]
red2blue = ["#FF0000", "#F90005", "#F4000A", "#EF000F", "#EA0014", "#E4001A", "#DF001F", "#DA0024", "#D50029",
            "#D0002E", "#CA0034", "#C50039", "#C0003E", "#BB0043", "#B60048", "#B0004E", "#AB0053", "#A60058",
            "#A1005D", "#9C0062", "#960068", "#91006D", "#8C0072", "#870077", "#82007C", "#7C0082", "#770087",
            "#72008C", "#6D0091", "#680096", "#62009C", "#5D00A1", "#5800A6", "#5300AB", "#4E00B0", "#4800B6",
            "#4300BB", "#3E00C0", "#3900C5", "#3400CA", "#2E00D0", "#2900D5", "#2400DA", "#1F00DF", "#1A00E4",
            "#1400EA", "#0F00EF", "#0A00F4", "#0500F9", "#0000FF"]
green2red = ["#00FF00", "#05F900", "#0AF400", "#0FEF00", "#14EA00", "#1AE400", "#1FDF00", "#24DA00", "#29D500",
             "#2ED000", "#34CA00", "#39C500", "#3EC000", "#43BB00", "#48B600", "#4EB000", "#53AB00", "#58A600",
             "#5DA100", "#629C00", "#689600", "#6D9100", "#728C00", "#778700", "#7C8200", "#827C00", "#877700",
             "#8C7200", "#916D00", "#966800", "#9C6200", "#A15D00", "#A65800", "#AB5300", "#B04E00", "#B64800",
             "#BB4300", "#C03E00", "#C53900", "#CA3400", "#D02E00", "#D52900", "#DA2400", "#DF1F00", "#E41A00",
             "#EA1400", "#EF0F00", "#F40A00", "#F90500", "#FF0000"]


def evenly_select(N, M):
    """
    Evenly select M elements from a list of length N.
    This returns a list of True/False (actually 0,1)
    See https://stackoverflow.com/questions/46494029/nearly-evenly-select-items-from-a-list

    Then use itertools.compress to create the new list

    :param N: The length of the list
    :param M: The number of elements to return
    :return: A list of 0/1 where 1 should be selected
    """
    if N == M:
        return np.ones(N, dtype=int)
    assert N > M
    if M > N/2:
        cut = np.ones(N, dtype=int)
        q, r = divmod(N, N-M)
        indices = [q*i + min(i, r) for i in range(N-M)]
        cut[indices] = False
    else:
        cut = np.zeros(N, dtype=int)
        q, r = divmod(N, M)
        indices = [q*i + min(i, r) for i in range(M)]
        cut[indices] = True

    return cut

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

def plotmap(ll, dd, outputfile, alpha, linewidth=1, bounds=None, maxdist=1, maxlinewidth=6, colorcontinents=False, plotintensity=False):
    """
    Plot the map of the dna distances and lat longs
    :param ll: The lon-lats
    :param dd: The distances to use
    :param outputfile: The file name to write the image to
    :param bounds: the boundary for the lat long as a 4-ple array
    :param maxdist: The maximum distance that we will scale to be maxlinewidth
    :param colorcontinents: color lines that go to different continents a different color (currently yellow)
    :param plotintensity: plot the colors of the lines by intensity vs. setting each number a color
    :return:
    """
    global verbose

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


    # note that now we calculate where everything should be and then plot it based on maximum values!
    dotat = {}
    for lonlat in ll.values():
        if bounds and ((lonlat[1] < bounds[0] or lonlat[1] > bounds[2]) or (lonlat[0] < bounds[1] or lonlat[0] > bounds[3])):
            if verbose:
                sys.stderr.write("Not in bounding box: {}\n".format(lonlat))
            continue
        if verbose:
            sys.stderr.write("Kept location: {}\n".format(lonlat))
        # plt.plot(lonlat[0], lonlat[1], 'o', color='Black', alpha=alpha, markersize=10, transform=ccrs.PlateCarree())
        dotat[(lonlat[0], lonlat[1])] = dotat.get((lonlat[0], lonlat[1]), 0) + 1

    maxdot = max(dotat.values())
    sys.stderr.write(f"Maximum dot density is {maxdot}\n")
    # we make the mean 50% intensity this time
    meandot = np.mean(list(dotat.values()))
    sys.stderr.write(f"The mean dot density is {meandot}\n")
    print()
    # now we color the dots based on the intensity of each dot!
    if 0:
        for tple in dotat:
            dotalpha = dotat[tple] / maxdot
            plt.plot(tple[0], tple[1], 'o', color='Black', alpha=dotalpha, markersize=10, transform=ccrs.PlateCarree())

    for tple in dotat:
        dotalpha = (dotat[tple] / meandot) * 0.5
        if dotalpha > 1:
            dotalpha = 1
        plt.plot(tple[0], tple[1], 'o', color='Black', alpha=dotalpha, markersize=10, transform=ccrs.PlateCarree())

    # how many lines and circles do we draw?
    circleat = {}
    circledata = {}
    lineat   = {}
    linedata = {}
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
            samecontinent = True
            if colorcontinents:
                # figure out if they are from the same continent
                m = re.search('\d{8}_(\w+)\_\d', idx1)
                cont1 = country2continent.get(m.groups(0)[0], "unknown")
                m = re.search('\d{8}_(\w+)\_\d', idx2)
                cont2 = country2continent.get(m.groups(0)[0], "unknown")
                if cont1 != cont2:
                    linecolor = 'yellow'
                    scaledalpha = alpha * 0.25
                    samecontinent = False

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

                #circ = Circle((circlon, circlat), transform=ccrs.Geodetic(), radius=radius,
                #                     linewidth=linewidth, alpha=scaledalpha, color=linecolor, fill=False)
                # ax.add_artist(circ)

                circleat[(circlon, circlat)] = circleat.get((circlon, circlat), 0) + 1
                circledata[(circlon, circlat)] = {
                    'radius' : radius,
                    'linewidth' : linewidth,
                    'alpha' : scaledalpha,
                    'color' : linecolor,
                    'fill'  : False
                }



            else:
                # plot a red line between two points
                #plt.plot([ll[idx1][0], ll[idx2][0]], [ll[idx1][1], ll[idx2][1]], color=linecolor, linewidth=linewidth,
                #         alpha=scaledalpha, transform=ccrs.Geodetic())

                linecoords = "\t".join(map(str, [ll[idx1][0], ll[idx2][0], ll[idx1][1], ll[idx2][1]]))

                lineat[linecoords] = lineat.get(linecoords, 0) + 1

                linedata[linecoords] = {
                    'x' : [ll[idx1][0], ll[idx2][0]],
                    'y' : [ll[idx1][1], ll[idx2][1]],
                    'color' : linecolor,
                    'linewidth' : linewidth,
                    'alpha' : scaledalpha,
                    'samecontinent' : samecontinent
                }


    # plot the circles
    circleval = sorted(set(circleat.values()))
    maxcircle = max(circleval)
    sys.stderr.write(f"The maximum circle is {maxcircle}\n")
    sys.stderr.write(f"There are {len(circleval)} circle values\n")
    # evenly select these colors from the list
    circlecolors = list(compress(red2blue, evenly_select(len(red2blue), len(circleval))))

    # do we want to do this by intensity or by number
    for tple in circleat:
        if plotintensity:
            idx = int((circleat[tple] / maxcircle) * len(red2blue)) - 1
            circlecolor = red2blue[idx]
        else:
            idx = circleval.index(circleat[tple])
            circlecolor = circlecolors[idx]

        circ = Circle((tple[0], tple[1]), transform=ccrs.Geodetic(), radius=circledata[tple]['radius'],
                             linewidth=circledata[tple]['linewidth'], alpha=circledata[tple]['alpha'], color=circlecolor, fill=circledata[tple]['fill'])
        ax.add_artist(circ)

    # plot the lines
    lineval = sorted(set(lineat.values()))
    maxline = max(lineval)
    sys.stderr.write(f"The maximum line is {maxline}\n")
    sys.stderr.write(f"There are {len(lineval)} line values\n")
    # evenly select these colors from the list
    linecolors = list(compress(red2blue, evenly_select(len(red2blue), len(lineval))))
    green2yellowline = list(compress(green2yellow, evenly_select(len(green2yellow), len(lineval))))
    for tple in lineat:
        if plotintensity:
            idx = int((lineat[tple] / maxline) * len(red2blue)) - 1
            if linedata[tple]['samecontinent']:
                colorline = red2blue[idx]
            else:
                colorline = green2yellow[idx]
        else:
            idx = lineval.index(lineat[tple])
            if linedata[tple]['samecontinent']:
                colorline = linecolors[idx]
            else:
                colorline = green2yellowline[idx]

        plt.plot(linedata[tple]['x'], linedata[tple]['y'], color=colorline, linewidth=linedata[tple]['linewidth'],
                 alpha=linedata[tple]['alpha'], transform=ccrs.Geodetic())



    #    plt.colorbar(CS3)

    # add a color bar for the lines and circles
    rect = Rectangle((10, 10), 140, 130, linewidth=5, edgecolor='b', facecolor='none')
    #plt.legend(handles=[rect])

    #grad = plt.imshow([[0.,1.], [0.,1.]], cmap = plt.cm.Reds, interpolation = 'bicubic')
    #plt.legend(handles=[grad])
    #fig = plt.figure()
    #ax2 = fig.add_axes()
    #ax2.add_patch(rect)


    #plt.show()
    plt.savefig(outputfile)

    # sys.stderr.write("We drew a max of {} circles\n".format(max(circleat.values())))
    # sys.stderr.write("And we drew a max of {} lines\n".format(max(lineat.values())))
    sys.stderr.write("Circles,{}\nLines,{}\n".format(",".join(map(str, circleat.values())), ",".join(map(str, lineat.values()))))
    sys.stderr.write("Dots,{}\n".format(",".join(map(str, dotat.values()))))

    sys.stderr.write("\nMAXIMUM VALUES\nCirlces: {}\nLines: {}\nDots: {}\n".format(max(circleat.values()),
                                                                                   max(lineat.values()),
                                                                                   max(dotat.values())
                                                                                   ))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot a map using ete and lat lon')
    parser.add_argument('-i', help='id.map file with lat/lon information', required=True)
    parser.add_argument('-m', help='cophenetic map file with same ids as id.map', required=True)
    parser.add_argument('-o', help='output file name', required=True)
    parser.add_argument('-a', help='alpha level for lines. Default=0.25', type=float, default=0.25)
    parser.add_argument('-l', help='linewidth for the lines connecting similar sites', default=1, type=float)
    parser.add_argument('-c', help='color the lines between continents yellow', action='store_true')
    parser.add_argument('-n', help='Plot the intensity as a fraction of max value', action='store_true')
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
    plotmap(lonlat, dist, args.o, args.a, linewidth=args.l, bounds=bounds, colorcontinents=args.c, plotintensity=args.n)
