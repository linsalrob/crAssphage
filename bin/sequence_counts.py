"""
Count information from a fasta file for Supplementary Table 4
"""

import os
import sys
import argparse


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
                seqid = line
                seq = ""
            else:
                seq += line
        seqs[seqid] = seq.upper()
    return seqs



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Count the sequences for supplementary table 4")
    parser.add_argument('-f', help='fasta file to count')
    args = parser.parse_args()

    seqs_by_id = read_fasta(args.f)
    seqs_by_seqs = set()
    seqs_by_seqs.update(seqs_by_id.values())
    locality = 0
    for l in seqs_by_id:
        if 'locality' in l:
            locality += 1

    print("{} / {} / {}".format(len(seqs_by_id), len(seqs_by_seqs), locality))