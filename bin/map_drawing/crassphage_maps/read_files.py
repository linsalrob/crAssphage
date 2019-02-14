import sys
import re
import gzip


def get_lon_lat(idf, maxtoget=50000, verbose=False):
    """
    Get the longitude and latitude of different ids. Note that we have longitude first to work with cartopy
    :param idf: the id.map file
    :param maxtoget: the maxiumum number of ids to get. This is just for debugging
    :param verbose: more output
    :return:
    """
    lonlat = {}
    count = 0
    with open(idf, 'r') as fin:
        for l in fin:
            if count > maxtoget:
                if verbose:
                    sys.stderr.write(f"Stopping because more than {maxtoget}\n")
                break
            count += 1
            s = re.search(r'latitude=(\S+)]', l)
            if not s:
                sys.stderr.write("No latitude in {}".format(l))
                continue
            lat = s.group(1)

            s = re.search(r'longitude=(\S+)]', l)
            if not s:
                sys.stderr.write("No longitude in {}".format(l))
                continue
            lon = s.group(1)
            p = l.split("\t")

            try:
                lat = float(lat)
                lon = float(lon)
            except:
                sys.stderr.write("There was an error parsing the latitude and longitude from {}\n".format(l))
                continue

            # lonlat[p[0]] = (lon, lat)
            newname = p[0].replace('|', '_')
            lonlat[newname] = (lon, lat)
    return lonlat


def closest_dna_dist(matrixfile, verbose=False):
    """
    Read the matrix file and get the id of the point with the closest distance that is not ourself
    :param matrixfile: The cophenetic matrix file to read
    :param verbose: more output
    :return: a dict of a node and its closest leaf
    """

    if verbose:
        sys.stderr.write("Getting closest distances\n")
    distances = {}

    if matrixfile.endswith('.gz'):
        with gzip.open(matrixfile, 'rt') as f:
            l = f.readline()
            ids = l.rstrip().split("\t")
            for i, name in enumerate(ids):
                if i == 0:
                    continue
                distances[name] = {}
            for l in f:
                data = l.rstrip().split("\t")
                for i, dist in enumerate(data):
                    if i == 0:
                        continue
                    distances[data[0]][ids[i]] = float(dist)
                    distances[ids[i]][data[0]] = float(dist)
    else:
        with open(matrixfile, 'r') as f:
            l = f.readline()
            ids = l.rstrip().split("\t")
            for i, name in enumerate(ids):
                if i == 0:
                    continue
                distances[name] = {}
            for l in f:
                data = l.rstrip().split("\t")
                for i, dist in enumerate(data):
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
