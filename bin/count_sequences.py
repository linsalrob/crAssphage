"""
Count the abundance of sequences in a fasta file
"""

import os
import sys
import argparse
from collections import Counter
from roblib import read_fasta

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Count the abundance of sequences in a fasta file")
    parser.add_argument('-f', help='fasta file', required=True)
    parser.add_argument('-v', help='verbose output', action="store_true")
    args = parser.parse_args()

    ids = read_fasta(args.f)
    lowercaseids = {k: v.lower() for k, v in ids.items()}
    seqs = Counter(lowercaseids.values())
    counts = Counter(seqs.values())
    for k,v in sorted(counts.items()):
        print(f"{k}\t{v}")


