"""
Rename our fasta files using locality:country:date:number format

Locality is probably (but not always) city, and matches the regexp [a-zA-Z0-9_] (i.e. \w+)
Country is the country and also matches the \w regexp
data is always an eight digit number e.g. 20171028 for October 28 2017 or 00000000 if we don't have the date
number is a number to disambiguate different samples from the same city on the same day.

The identifier can be separated with either:
    split(":")
    or
    \w+:\w+:\d{8}:\d+

"""

import os
import sys
import argparse
import re

def stripout(txt):
    """
    Replace white space and non \w characters in text
    :param txt: some text to clean
    :return: clean text
    """
    txt = txt.replace(' ', '_')
    return re.sub('\W', '', txt)


def rename(fafile, idmapfile, outputfa):
    """
    Rename the sequences
    :param fafile: input fasta file
    :param idmapfile: output id.map file with new name and original name
    :param outputfa: output fasta file
    :return: None
    """

    counts = {}

    idout = open(idmapfile, 'w')
    faout = open(outputfa, 'w')
    fatal_error = False
    notfirstsequence = False # we don't want to print a new line for the first sequence
    with open(fafile, 'r') as fa:
        for l in fa:
            l=l.strip()
            lcl = 'NA'
            cnt = 'NA'
            dt = '00000000'
            if l.startswith('>'):
                if notfirstsequence:
                    faout.write("\n")
                else:
                    notfirstsequence = True
                l=l.replace('>', '', 1)
                if 'locality' in l:
                    m = re.search('\[locality=(.*?)\]', l)
                    if not m:
                        sys.stderr.write("Even though its there, can't parse locality in {}\n".format(l))
                        fatal_error = True
                    lcl = stripout(m.groups()[0])
                else:
                    sys.stderr.write("FATAL ERROR: Locality was not found in {}\n".format(l))
                    fatal_error = True

                if 'country' in l:
                    m = re.search('\[country=(.*?)\]', l)
                    if not m:
                        sys.stderr.write("Even though its there, can't parse country in {}\n".format(l))
                        fatal_error = True
                    cnt = stripout(m.groups()[0])
                else:
                    sys.stderr.write("FATAL ERRROR: Country was not found in {}\n".format(l))
                    fatal_error = True


                if 'date' in l:
                    m = re.search('\[date=(.*?)\]', l)
                    if not m:
                        sys.stderr.write("Even though its there, can't parse country in {}\n".format(l))
                        fatal_error = True
                    dt = m.groups()[0]
                    if re.search('\D', dt):
                        sys.stderr.write("WARNING: the date in {} has non-digit entries\n".format(l))
                        dt = re.sub('\D', '', dt)
                    if not re.match('\d{8}', dt):
                        sys.stderr.write("WARNING the date in {} is not correct\n".format(l))
                        while not re.match('\d{8}', dt):
                            dt += '0'

                identifier = "{}|{}|{}|".format(lcl, dt, cnt)

                counts[identifier] = counts.get(identifier, 0) + 1
                identifier += str(counts[identifier])

                idout.write("{}\t{}\n".format(identifier, l))
                faout.write(">{}\n".format(identifier))
            else:
                faout.write(l)
    faout.write("\n")

    idout.close()
    faout.close()

    if fatal_error:
        sys.exit(-1)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-f', help="fasta file to rename", required=True)
    parser.add_argument('-o', help='output fasta file', required=True)
    parser.add_argument('-i', help='output id.map file that has the original information and the new information', required=True)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    rename(args.f, args.i, args.o)
