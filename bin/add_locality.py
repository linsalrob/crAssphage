"""

Read a fasta file and add a locality location to every entry that does not have it.

"""

import os
import sys
import argparse
import re
import geocoding


def parse_fasta(fafile, outfile):
    """
    Parse the fasta file
    """

    changed=False
    with open(outfile, 'w', encoding='utf-8') as out:
        with open(fafile, 'r') as f:
            for l in f:
                if not l.startswith('>'):
                    out.write(l)
                    continue
                if 'latlon' not in l:
                    sys.stderr.write("ERROR: No latlon in {}".format(l))
                    out.write(l)
                    continue

                m = re.search('latlon=(.*?)\]', l)
                lat, lon = m.groups()[0].split(',')
                locality, country = geocoding.latlon_to_place(lat, lon)

                newline = l.rstrip()
                if locality and 'locality' not in l:
                    newline = newline + " [locality={}]".format(locality)
                    changed = True
                if country and 'country' not in l:
                    newline = newline + " [country={}]".format(country)
                    changed = True
                out.write("{}\n".format(newline))
    
    if changed:
        sys.stderr.write("Changes were made to {}\n".format(fafile))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Use geocode to add a locality to the metadata')
    parser.add_argument('-f', help='fasta file to add locality to', required=True)
    parser.add_argument('-o', help='fasta output file to write', required=True)
    args = parser.parse_args()

    parse_fasta(args.f, args.o)






