#__python__

import os
import sys

try:
    ccf = sys.argv[1]
    gcf = sys.argv[2]
except:
    sys.exit(sys.argv[0] + " <country counts file> <countries list>")


if not os.path.exists(ccf):
    sys.exit(ccf + " not found")
if not os.path.exists(gcf):
    sys.exit(gcf + " not found")

country={}
with open(gcf, 'r') as gc:
    for l in gc:
        if l.startswith('#'):
            continue
        p=l.strip().split("\t")
        if len(p) < 3:
            continue
        country[p[3]] = [p[1], p[2], p[0]]

with open(ccf, 'r') as cc:
    for l in cc:
        if l.startswith('#'):
            continue
        p=l.strip().split("\t")
        if p[0] not in country:
            sys.stderr.write("Country: {} not found\n".format(p[0]))
            for c in country:
                if c.startswith(p[0]):
                    sys.stderr.write("\tPerhaps: {}\n".format(c))




