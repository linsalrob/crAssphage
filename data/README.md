# Data

This directory contains data that you might need to analyze the crAssphage sequences, including the original submission to GenBank that was subsequently deleted. 

The JQ995537 files are the sequence files for crAssphage. JQ995537.fna is the nucleic acid fasta sequence of the crAssphage genome. The other three files are the blast formated databases for that genome.

[countries.csv](countries.csv) is a list of country codes and lat/lon locations. This is taken from the [Google Developers Toolkit](https://developers.google.com/public-data/docs/canonical/countries_csv) and you should use that version as it is more up to date. I added a couple of redundancies for additional spelling.

[localities.tsv](localities.tsv) is a tab separated list of lat/lon locations that we have samples for and the locality and country that sample comes from. We have the locality and country data in two formats: the original UTF-8 UniCode encoding, and a regular ASCII encoding that we use in the sequence files (as many bioinformatics programs still can't handle UTF-8).

[localities.db](localities.db) is an SQLite3 database of the same data. Note that you should check the respective update times for these two files, as we have code to automatically update the database as we add new locations, but that doesn't necessarily propagate into the tab separated text file above. This database currently has a single table, location, with the attributes [*latitude*, *longitude*, *locality*, *country*, *ascii_locality*, *ascii_country*]. 

Note that the locality data is taken from the Google Maps API but then manually curated for missing elements and to create the ASCII versions of place names. Any errors are because of Rob!
