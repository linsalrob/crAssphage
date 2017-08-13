import argparse
import os
import sys
import re
from roblib import sequences

"""

Add information about the country, date, and lat lon to the metagenome reads.

In the metadata file, country is column 51, collection date is column 29, and 
lat lon is column 80.

"""

def parse_sample_ids(sampleidfile):
    """
    Read the sample ID file and return a hash
    of SRS to list of SRR (SRS->[SRR])
    """

    srs={}
    with open(sampleidfile, 'r') as f:
        for l in f:
            p=l.strip().split("\t")
            if p[1] not in srs:
                srs[p[1]]=[]
            srs[p[1]].append(p[0])

    return srs


def parse_metadata(metadatafile, srs):
    """
    Given the metadata file name and the srs hash (from parse sample IDs)
    we construct a read hash that has every read and country, date, latlon.

    Empty values are labeled as 'unknown'
    """

    srr={}
    with open(metadatafile, 'r') as f:
        for l in f:
            if l.startswith('Sample Accession'):
                continue
            p=l.split("\t")
            m=re.search('SRA\|(\w+);', p[1])
            if not m:
                sys.stderr.write("Couldn't find an SRA sample in {}".format(l))
                continue
            srsid = m.groups()[0]
            if srsid not in srs:
                sys.stderr.write("No SRS info found for {}\n".format(srsid))
                continue
            if ':' in p[51]:
                countryparts = p[51].split(':')
                p[51] = countryparts[0]
            if p[51] and p[29]:
                for srrid in srs[srsid]:
                    srr[srrid]=[p[51], p[29], p[80]]

    return srr

def parse_fasta(fastafile, srr):
    """
    Given the srr hash and the fasta file add the data from one to the other!
    """

    for (seqid, seq) in sequences.stream_fasta(fastafile):
         m=re.match('([^_]*)', seqid)
         srrid=m.groups()[0]
         if srrid in srr:
             sys.stdout.write(">{} [country={}] [date={}] [lat_lon={}]\n{}\n".format(
                 seqid, srr[srrid][0], srr[srrid][1], srr[srrid][2], seq))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Find just those seqs with a country, date, latlon and print them')
    parser.add_argument('-f', help='fasta file', required=True)
    parser.add_argument('-s', help='file with list of sample ids and reads (sample_ids.txt)', required=True)
    parser.add_argument('-m', help='file with country information (crAssphage.country.date.tsv)', required=True)
    args = parser.parse_args()

    srs = parse_sample_ids(args.s)
    srr = parse_metadata(args.m, srs)
    parse_fasta(args.f, srr)
            


