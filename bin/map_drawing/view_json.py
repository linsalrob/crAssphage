"""
view the json file data structure
"""

import os
import sys
import argparse
import json
import pprint

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="print the contents of a json file")
    parser.add_argument('-f', help='Json file')
    parser.add_argument('-v', help='verbose output', action="store_true")
    args = parser.parse_args()

    p = pprint.PrettyPrinter(indent=4)
    j = json.load(open(args.f, 'r'))
    p.pprint(j)
