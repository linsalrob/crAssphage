"""
Create a tab separated file of metadata with column headers that start
with a #.

We read the id map file and make an entry for every sequence, making sure
to check for missing data
"""

import os
import sys
from collections import OrderedDict
import argparse
import re

__author__ = 'Rob Edwards'



def parse_fasta(idfile):
    """
    Parse the fasta file and return a hash of sequence ids and metadata
    
    We also do some sanity checks on the data
    """
    seqids = set()
    data = OrderedDict()
    fatalerror = False

    with open(idfile, 'r') as f:
        for l in f:
            p=l.strip().split("\t")
            seqids.add(p[0])
            for meta in re.findall('\\[(.*?)\\]', p[1]):
                tag,value=meta.split('=')
                if tag not in data:
                    data[tag] = {}
                data[tag][p[0]]=value


        # now we have read the file, lets check a few things out
        
        wantedkeys = ['latitude', 'longitude', 'latlon', 'locality', 'country']
        for w in wantedkeys:
            if w not in data:
                sys.stderr.write("FILE ERROR: No {} entries were found in {}\n".format(w, idfile))
                continue
            for s in seqids:
                if s not in data[w]:
                    sys.stderr.write("ENTRY ERROR: No {} was found for {}\n".format(w, s))

    return seqids, data


def print_table(seqids, data, outputfile, separator='\t'):
    """
    Write the metadata to a tab separated text file
    """

    tags = data.keys()
    with open(outputfile, 'w') as out:
        out.write(separator.join(["#Sequence ID"] + list(tags)) + "\n")
        for s in seqids:
            out.write(s)
            for t in tags:
                out.write("{}{}".format(separator, data[t].get(s, "")))
            out.write("\n")






if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Make a metadata file from the id.map file')
    parser.add_argument('-i', help='id.map file to make the metadata from', required=True)
    parser.add_argument('-o', help='output file', required=True)
    args = parser.parse_args()

    seqids, data = parse_fasta(args.i)
    print_table(seqids, data, args.o)

