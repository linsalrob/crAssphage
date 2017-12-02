"""
Read a cophenetic matrix
"""

import os
import sys
import argparse
import gzip

def pairwise_distances(matrixfile):
    """
    Read the matrix file and return a hash of all distances
    :param treefile: The cophenetic matrix file to read
    :return: a dict of each entry and its distances
    """

    global verbose

    distances = {}

    f = None
    if matrixfile.endswith('.gz'):
        f = gzip.open(matrixfile, 'rt')
    else:
        f = open(matrixfile, 'r')

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

    f.close()

    return distances

def closest_dna_dist(matrixfile):
    """
    Read the matrix file and get the id of the point with the closest distance that is not ourself
    :param treefile: The cophenetic matrix file to read
    :return: a dict of a node and its closest leaf
    """

    global verbose
    distances = {}

    f = None
    if matrixfile.endswith('.gz'):
        f = gzip.open(matrixfile, 'rt')
    else:
        f = open(matrixfile, 'r')

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

    return closest

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('-f', help='')
    parser.add_argument('-v', help='verbose output')
    args = parser.parse_args()