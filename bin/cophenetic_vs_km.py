"""
Generate all pairwise comparisons of cophenetic distance and physical distance
"""

import os
import sys
import argparse
import math

from scipy.stats.stats import pearsonr
import matplotlib.pyplot as plt

from km_distance import get_lon_lat, latlon2distance
from cophenetic import pairwise_distances, closest_dna_dist


def write_distances(dist, lonlat, outputfile, writeall=False):
    """
    Write the distances to a file
    :param dist: the distance hash
    :param lonlat: the latitude and longitude hash
    :param outputfile: the file to write
    :param writeall: whether to write them all or just the unique set
    :return: nothing
    """
    uniqdist = {}
    with open(outputfile, 'w') as out:
        out.write("km\tGenetic Distance\n")
        for i in dist:
            for j in dist[i]:
                ii = i.replace('|', '_')
                jj = j.replace('|', '_')
                km = latlon2distance(lonlat[ii][0], lonlat[ii][1], lonlat[jj][0], lonlat[jj][1])
                if writeall:
                    out.write("{}\t{}\t{}\t{}\n".format(i, j, km, dist[i][j]))
                else:
                    if km not in uniqdist:
                        uniqdist[km] = set()
                    uniqdist[km].update([dist[i][j]])
        if not writeall:
            for k in uniqdist:
                for d in uniqdist[k]:
                    out.write("{}\t{}\n".format(k,d))

def distances_to_list(dist, latlon):
    """
    Create two arrays, one of the km, and one of the genetic distance
    :param dist: the distance hash
    :param latlon: the lat lon has
    :return: [km array] [genetic distance array]
    """
    km = []
    gd = []

    for i in dist:
        for j in dist[i]:
            ii = i.replace('|', '_')
            jj = j.replace('|', '_')
            km.append(latlon2distance(lonlat[ii][0], lonlat[ii][1], lonlat[jj][0], lonlat[jj][1]))
            gd.append(dist[i][j])

    return km, gd

def pearson_correlation(km, gd):
    """
    Calculate the pearson correlation between the km and gd
    :param km: the km array
    :param gd: the gd array
    :return: the pearson correlation and p value
    """

    return pearsonr(km, gd)


def plot_km_gd(km, gd, outputfile):
    """
    Plote the figure in outputfile
    :param km: the km array
    :param gd: the gd array
    :param outputfile: the picture file to write
    :return:
    """

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(km, gd, 'kx')

    ax.set_xlabel("Physical Distance (km)")
    ax.set_ylabel("Genetic Distance")

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()


    fig.set_facecolor('white')

    plt.tight_layout()
    plt.savefig(outputfile)

def settoepsilon(arr, e=1e-10):
    arr = [e if 0 == x else x for x in arr]
    return arr



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot a map using ete and lat lon')
    parser.add_argument('-i', help='id.map file with lat/lon information', required=True)
    parser.add_argument('-m', help='cophenetic map file with same ids as id.map', required=True)
    parser.add_argument('-o', help='file name to write the matrices to')
    parser.add_argument('-a', help='output all pairwise distances. Otherwise just a unique set', action='store_true')
    parser.add_argument('-p', help='file name for the plot of km vs genetic distance')
    parser.add_argument('-c', help="Calcuate the pearson correlation", action='store_true')
    parser.add_argument('-l', help='Plot log of genetic distance', action='store_true')
    parser.add_argument('-x', help='maximum genetic distance to consider', type=float)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    global verbose
    verbose = False
    if args.v:
        verbose = True

    lonlat = get_lon_lat(args.i)
    dist = pairwise_distances(args.m)
    #dist = closest_dna_dist(args.m)

    if args.x:
        tempdist = {x:{} for x in dist}
        for x in dist:
            for y in dist[x]:
                if dist[x][y] < args.x:
                    tempdist[x][y] = dist[x][y]
        dist=tempdist



    if args.o:
        write_distances(dist, lonlat, args.o, args.a)

    if not args.p and not args.c:
        sys.exit(0)

    km, gd = distances_to_list(dist, lonlat)

    if args.l:
        gd = settoepsilon(gd, 1e-5)
        gd = [math.log(x)/math.log(10) for x in gd]
        km = settoepsilon(km, 1)
        km = [math.log(x)/math.log(10) for x in km]

    if args.c:
        pe, pv = pearsonr(km, gd)
        print("Pearson correlation: {} p-value: {}".format(pe, pv))

    if args.p:
        plot_km_gd(km, gd, args.p)


