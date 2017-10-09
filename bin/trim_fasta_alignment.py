"""
Trim a fasta alignment file based on the columns with - and then remove sequences that don't have enough coverage.
We assume all the sequences are the same length
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
        seqs[seqid] = seq
    return seqs




if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Trim fasta format alignment files based on gaps and lengths')
    parser.add_argument('-f', help='fasta file of alignment', required=True)
    parser.add_argument('-c', help='minimum coverage of a COLUMN. Default = 0.9. Columns < than this will be removed', default=0.9, type=float)
    parser.add_argument('-r', help='minimum coverage of a ROW (sequence). Default = 0.9. Sequences < than this will be removed', default=0.9, type=float)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()


    coltokeep = [] # which columns are we going to output
    rowtokeep = [] # which rows do we output


    dna = read_fasta(args.f)
    ids = list(dna.keys())
    numseqs = len(ids)

    for i in range(len(dna[ids[0]])):
        count = 0
        for seq in ids:
            if '-' != dna[seq][i]:
                count += 1
        frac = 1.0 * count/numseqs
        if 1.0 * count/numseqs > args.c:
            coltokeep.append(i)

    if len(coltokeep) == 0:
        sys.exit("We did not elect to keep any bases in the alignment!")

    finalbases = len(coltokeep)
    for seq in ids:
        count = 0
        for i in coltokeep:
            if "-" != dna[seq][i]:
                count += 1
        if 1.0 * count/finalbases > args.r:
            rowtokeep.append(seq)
            sys.stderr.write("Accepted sequence {} because it has {}% bases [count: {} finalbases: {}]\n".format(seq, 1.0 * count / finalbases, count, finalbases))
        elif args.v:
            sys.stderr.write("Rejected sequence {} because it has {}% bases\n".format(seq, 1.0 * count/finalbases))

    for seq in rowtokeep:
        print(seq)
        for i in coltokeep:
            sys.stdout.write(dna[seq][i])
        sys.stdout.write("\n")
