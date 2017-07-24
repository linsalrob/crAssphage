# Code for analysing the data.

This is some code that we use frequently in the analysis of the sequences, and is included here for convenience. 

All of this code is released as-is, and under the [MIT Licence](../../LICENSE). You are welcome to use it or ignore it!


## Code for making the trees

The following code is used as part of the Makefile for building trees

* `renumber_fasta.pl`
This takes a fasta file and changes the IDs to seequential numbers and makes a file called id.map that has the new and old IDs.
* `trim_alignment.py`
This code removes "columns" in the alignment if they have too many gaps
* `negative_branch_lengths.py`
This code adjusts negative branch lengths, setting them to zero and moving the difference from zero to the adjacent branch as [discussed here](http://www.icp.ucl.ac.be/~opperd/private/neighbor.html)
* `newick.py`
This is a python module for parsing newick trees used for renaming leaves and adjusting negative branch lengths.
* `phylip2clustal.py`
This is a simple biopython module for converting phylip format trees (Which we use as the default, since we use phylip for the tree) to clustal format.
* `rename_trees.py`
This code changes the names on the leaves of the trees. 


## Other code

* `idmap2distance.py`
This code will read the id.map file and calculate the distance between every sampled site. It generates a matrix of physical distance in km between every site.
* `dnadist2anova.py`
This code reads the seqs.dnadist file and generates a matrix of dna distances for all the samples. It also adds secondary information to those samples.

