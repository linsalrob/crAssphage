"""
A revised version of rename trees that uses the ete3 library to parse the tree
"""

import os
import sys
import argparse
import re
from ete3 import Tree


if __name__ == '__main__':

    tags = ['address', 'altitude', 'country', 'date', 'id', 'latitude', 'location', 'longitude', 'name', 'note', 'samplecount', 'site', 'volunteer']

    parser = argparse.ArgumentParser(description='Parse a tree')
    parser.add_argument('-t', help='tree file', required=True)
    parser.add_argument('-i', help='id map file', required=True)
    parser.add_argument('-n', help='name of the tag to use in the label. You can provide multiple tags and they will be joined with _ Valid tags are {}'.format(tags), action='append')
    args = parser.parse_args()

    if not args.n:
        args.n = ['id']

    for n in args.n:
        if n.lower() not in tags:
            sys.exit('Sorry {} is not a valid tag. It must be one of {} '.format(args.n, tags))

    idmap = {}
    tagcount = {}
    with open(args.i, 'r') as f:
        for l in f:
            p=l.strip().split("\t")
            # note that some programs replace | with _ so we use both
            altid = p[0].replace('|', '_')
            tagarr = []
            if 'id' in args.n:
                tagarr.append(p[1].split()[0])
            else:
                m = re.findall('\[(\S+)\=(.*?)\]', p[1])
                tagdata = {t[0]:t[1] for t in m}
                for n in args.n:
                    if n.lower() in tagdata:
                        tagarr.append(tagdata[n.lower()])
            tagstring = "_".join(tagarr)

            if '' == tagstring:
                sys.stderr.write("No {} data found in sample: {}\n".format(args.n, l.strip()))
                idmap[p[0]] = p[1].split()[0]
                idmap[altid] = p[1].split()[0]
            else:
                tagcount[tagstring] = tagcount.get(tagstring, 0)+1
                idmap[p[0]] = "{}_{}".format(tagstring, tagcount[tagstring])
                idmap[altid] = "{}_{}".format(tagstring, tagcount[tagstring])

    tree = Tree(args.t)
    for leaf in tree:
        oldname = leaf.name
        oldname = oldname.replace('_R_', '')
        if oldname not in idmap:
            sys.stderr.write("ERROR: {} is not in {} as an ID\n".format(oldname, args.i))
            continue
        leaf.name = idmap[oldname]

    print(tree.write())
