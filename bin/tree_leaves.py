"""
Very simple, just report the number of leaves in a tree.
"""

import os
import sys
import argparse
from ete3 import Tree

def print_leaves(treefile):
    """
    Print how many leaves there are in this tree
    :param treefile: The tree file to parse
    :return:
    """

    count=0
    tree = Tree(treefile)
    for leaf in tree:
        count+=1
    print("{} : {} leaves".format(treefile, count))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Count tree leaves not tea leaves")
    parser.add_argument('-f', help='Tree to count leaves in')
    parser.add_argument('-v', help='verbose output', action="store_true")
    args = parser.parse_args()

    print_leaves(args.f)