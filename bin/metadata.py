"""
Create a tab separated file of metadata with column headers that start
with a #.

We read the fasta file and make the file for every sequence, making sure
to check for missing data
"""

import os
import sys
from collections import OrderedDict
import argparse
import re

__author__ = 'Rob Edwards'



def parse_fasta(fafile):
    """
    Parse the fasta file and return a hash of sequence ids and metadata
    
    We also do some sanity checks on the data
    """
    
    seqids = set()
    data = OrderedDict()
    fatalerror = False

    with open(fafile, 'r') as f:
        for l in f:
            if not l.startswith('>'):
                continue
            seqid=l.split(' ')[0].replace('>', '', 1)
            if (seqid in seqids):
                sys.stderr.write("FATAL: Duplicate sequence ID: {} in {}\n".format(seqid, fafile))
                fatalerror = True
            seqids.add(seqid)
            for meta in re.findall('\\[(.*?)\\]', l):
                tag,value=meta.split('=')
                if tag not in data:
                    data[tag] = {}
                data[tag][seqid]=value


        # now we have read the file, lets check a few things out
        
        wantedkeys = ['latitude', 'longitude', 'latlon', 'city', 'country']
        for w in wantedkeys:

            if w not in data:
                sys.stderr.write("FILE ERROR: No {} entries were found in {}\n".format(w, fafile))
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
    parser = argparse.ArgumentParser(description='Make a metadata file from a fasta file')
    parser.add_argument('-f', help='fasta file to make the metadata from', required=True)
    parser.add_argument('-o', help='output file', required=True)
    args = parser.parse_args()

    seqids, data = parse_fasta(args.f)
    print_table(seqids, data, args.o)

