"""
Tree collapser.

Collapse nodes in a phylogenetic tree.

(c) 2018 Alessandro Rossi
"""

import os
import sys
import argparse
from ete3 import Tree


def collapse_nodes(inputfile, threshold, output_name):
    """
    Collapse nodes less than threshold
    :param inputfile: The tree file
    :param threshold: The threshold on which to collapse
    :param output_name: The output file name
    :return:
    """
    input_tree = Tree(inputfile)

    for node in input_tree.get_descendants():
        if node.support < threshold and node.is_leaf() == False:
            node.delete(preserve_branch_length=True)

    Tree.write(input_tree, outfile=output_name)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Collapse tree nodes below a certain threshold")
    parser.add_argument('-f', help='Input tree file', required=True)
    parser.add_argument('-t', help='Threshold for collapsing the nodes (float)', required=True, type=float)
    parser.add_argument('-o', help='output tree file', required=True)
    parser.add_argument('-v', help='verbose output', action="store_true")
    args = parser.parse_args()

    collapse_nodes(args.f, args.t, args.o)