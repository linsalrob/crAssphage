"""
From the id.map file create a matrix of samples and their (physical)
distances based on the latitude/longitude. 

We convert lat/lon to a km distance and then make a pairwise table.
"""

import os
import sys
import argparse
import re
import math

__author__ = 'Rob Edwards'



def latlon2distance(lat1, long1, lat2, long2, miles=False):
    """Convert two coordinates to distance.

    This is an approximation since the earth is not spherical, but accuracy is <100m, especially for close points

    This code was taken from http://www.johndcook.com/python_longitude_latitude.html

    Latitude is measured in degrees north of the equator; southern locations have negative latitude.
    Similarly, longitude is measured in degrees east of the Prime Meridian. A location 10deg west of
    the Prime Meridian, for example, could be expressed as either 350deg  east or as -10deg east.

    Arguments: lat1, long1; lat2, long2; miles is a boolean. If you want miles set it to true. Else set it to false

    """

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
    arc = math.acos(cos)

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

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Make a (geographic) distance matrix from the id.map file')
    parser.add_argument('-i', help='id.map output file from the renumbering code', required=True)
    parser.add_argument('-o', help='output file', required=True)
    args = parser.parse_args()

    loc = {}
    with open(args.i, 'r') as f:
        for l in f:
            p = l.strip().split("\t")
            m=re.findall('^(\S+)', p[1])
            sampleid = m[0]
            m=re.search('\[latitude=(.*?)\]', p[1])
            lat=m.group(1)
            m=re.search('\[longitude=(.*?)\]', p[1])
            lon=m.group(1)
            loc[sampleid]=[float(lat), float(lon)]

    names = list(loc.keys())
    names.sort()

    with open(args.o, 'w') as out:
        out.write("sampleid\t" + "\t".join(names))
        out.write("\n")
        for i in range(len(names)):
            fn = names[i]
            out.write(fn)
            for j in range(len(loc)):
                tn = names[j]
                if fn == tn:
                    out.write("\t0")
                    continue
                dist = latlon2distance(loc[fn][0], loc[fn][1], loc[tn][0], loc[tn][1])
                out.write("\t{}".format(dist))
            out.write("\n")

