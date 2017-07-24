import os
import sys
import argparse
import os
import string
import sys

import re

__author__ = 'Rob Edwards'

"""
Python code to read a tree file, rename all the leaves based on a tab separated text file that you provide with the
-i option, and print out the new tree.


Unless you have a large tree, you should use rename_trees.py. That code parses the tree and takes each leaf and renames
it. In contrast, this version uses regular expressions to just identify the labels in the tree and rename them.  

"""


def clean_name(name):
    """
    Just clean out non-allowable characters in the name

    :param name: the new name
    :type name: str
    :return: the revised name
    :rtype: str
    """

    allowable = set(string.ascii_letters)
    allowable.update(set(string.digits))
    allowable.update({'_','-',':'})
    name.replace(' ', '_')
    return "".join(filter(lambda x: x in allowable, name))

def rename_leaf(name, idmap):
    """
    Rename the nodes of a tree based on id map

    :param root: the root node of the tree
    :type root: Node
    :param idmap: the id map
    :type idmap: dict
    :return: the new root node
    :rtype: Node
    """

    firstchar = name[0]
    name = name[1:]

    if name.startswith('_R_'):
        # this was reverse complemented by MAFFT
        name = name.replace('_R_', '')

    if name and name in idmap:
        name = clean_name(idmap[name])

    return firstchar + name


if __name__ == '__main__':

    tags = ['id', 'address','altitude','country','date','latitude','longitude','name','note', 'site']

    parser = argparse.ArgumentParser(description='Parse a tree')
    parser.add_argument('-t', help='tree file', required=True)
    parser.add_argument('-i', help='id map file', required=True)
    parser.add_argument('-n', help='name of the tag to use in the label. Valid tags are {}'.format(tags), default='id')
    args = parser.parse_args()

    if args.n.lower() not in tags:
        sys.exit('Sorry {} is not a valid tag. It must be one of {} '.format(args.n, tags))

    idmap = {}
    tagcount = {}
    with open(args.i, 'r') as f:
        for l in f:
            p=l.strip().split("\t")
            if tags == 'id':
                idmap[p[0]]=p[1].split()[0]
            else:
                m = re.findall('\[(\S+)\=(.*?)\]', p[1])
                tagdata = {t[0]:t[1] for t in m}
                if args.n.lower() in tagdata:
                    tagcount[tagdata[args.n]] = tagcount.get(tagdata[args.n], 0)+1
                    idmap[p[0]] = "_".join(map (str, [tagdata[args.n], tagcount[tagdata[args.n]]]))
                else:
                    if args.n.lower() != 'id':
                        sys.stderr.write("No {} data found in sample: {}\n".format(args.n, l.strip()))
                    idmap[p[0]] = p[1].split()[0]

    new = ""
    with open(args.t, 'r') as tre:
        for l in tre:
            new = new + re.sub('[(,]\w+', lambda m: rename_leaf(m.group(), idmap), l)
    print(new)

