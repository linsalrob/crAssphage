import os
import sys
import location_db as ldb

# create the database
cursor = ldb.create_database()

with open('../data/latlon_location.tsv', 'r', encoding='utf-8') as f:
    for l in f:
        p=l.strip().split("\t")
        ldb.save_location(p[0], p[1], p[2], p[3], p[4], p[5], cursor)



