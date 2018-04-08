"""
Create a tab separated file of metadata with column headers that start
with a #.

We read the fasta file and the id map file and make an entry for every sequence, making sure
to check for missing data


"""

import os
import sys
from collections import OrderedDict
import argparse
import re

__author__ = 'Rob Edwards'


def read_fasta(fname):
    """
    Read a fasta file and return a hash.
    :param fname: The file name to read
    :return: dict
    """

    seqs = {}
    seq = ""
    seqid = ""
    with open(fname, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if "" != seqid:
                    seqs[seqid] = seq
                seqid = line.replace('>', '', 1)
                seq = ""
            else:
                seq += line
        seqs[seqid] = seq
    return seqs




def parse_id_map(idmapfile, seqs, primer):
    """
    Parse the fasta file and return a hash of sequence ids and metadata
    
    We also do some sanity checks on the data
    :param idmapfile: the fasta file to check
    :param primer: the primer used for the PCR
    :return: dict, dict
    """

    data = OrderedDict()
    for x in ['seqid', 'primer', "name", "date", "yearmonth", "lat_lon", "locality", "latitude", "longitude", "latlon", "altitude", "address", "country", "location", "contact", "description", "method", "note", "bin_id", "sample", "sample_date", "sampleid", "sample_number", "sample_processing", "sampletype", "sex", "site", "source", "type", "university", "volunteer"]:
        data[x]={}

    fatalerror = False

    with open(idmapfile, 'r') as f:
        for l in f:
            seqid=l.strip().split("\t")[0]
            if seqid not in seqs:
                sys.stderr.write("Sequence id {} was not found in the sequence file\n".format(seqid))
                fatalerror = True

            data["seqid"][seqid]=seqid
            data['primer'][seqid]=primer
            for meta in re.findall('\[(.*?)\]', l):
                tag,value=meta.split('=')
                tag = tag.lower()
                if 'primer' == tag:
                    if primer != value:
                        sys.stderr.write("You told me the primer was {} but sample {} thinks the primer is {}\n".format(primer, l, value))
                    continue
                if 'notes' == tag:
                    tag = 'note'
                if tag not in data:
                    sys.stderr.write("Added a tag ({}) that was not in our defined list. Order may be different for primers A, B, and C\n".format(tag))
                    data[tag] = {}
                data[tag][seqid]=value


    # now we have read the file, lets check a few things out
    wantedkeys = ['latitude', 'longitude', 'latlon', 'locality', 'country']
    for w in wantedkeys:
        if w not in data:
            sys.stderr.write("FILE ERROR: No {} entries were found in {}\n".format(w, idmapfile))
            continue
        c=0
        for s in seqs:
            c+=1
            if c > 10:
                break
            if s not in data[w]:
                sys.stderr.write("ENTRY ERROR: No {} was found for {}\n".format(w, s))

    return data


def print_table(seqs, data, outputfile, separator='\t'):
    """
    Write the metadata to a tab separated text file
    """

    tags = [x for x in data.keys() if x != 'seqid']
    with open(outputfile, 'w') as out:
        out.write(separator.join(["#sequence ID", "alternate ID"] + tags + ["sequence\n"]))
        for s in seqs:
            tempid = data['seqid'][s].replace('|', '_')
            out.write("{}\t{}".format(tempid, data['seqid'][s]))
            for t in tags:
                out.write("{}{}".format(separator, data[t].get(s, "")))
            out.write("{}{}\n".format(separator, seqs[s].upper()))




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Make a metadata file that inlcudes all metadata and the DNA sequence from a fasta file')
    parser.add_argument('-f', help='fasta file to read the sequences from', required=True)
    parser.add_argument('-i', help='idmap to make the metadata from', required=True)
    parser.add_argument('-p', help='primer used for this data', required=True)
    parser.add_argument('-o', help='output file', required=True)
    args = parser.parse_args()

    seqs = read_fasta(args.f)
    data = parse_id_map(args.i, seqs, args.p)
    print_table(seqs, data, args.o)

