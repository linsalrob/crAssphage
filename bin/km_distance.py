"""
Read id maps and get distances
"""

import os
import sys
import argparse
import re
import math

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
            lonlat[p[0]] = (float(lon), float(lat))
            newname = p[0].replace('|', '_')
            lonlat[newname] = (float(lon), float(lat))

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



