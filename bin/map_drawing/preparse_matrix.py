"""
Parse the matrix file and write a JSON file for the data
"""

import os
import sys
import argparse
from roblib import bcolors
import json

from crassphage_maps import closest_dna_dist

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-m', help='matrix file of the cophenetic matrix to parse', required=True)
    parser.add_argument('-o', help='outputfile to write', required=True)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    dist = closest_dna_dist(args.m)
    with open(args.o, 'w') as jsonout:
        json.dump(dist, jsonout)
