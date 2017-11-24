# Trees and alignments of all the data.

The trees and alignments are made in several steps, although all the steps are listed in the [Makefile](PrimerA/Makefile)

First, we collect all the sequenes in fasta format. Hopefully, those sequences are symbolic links to the original data directories (they are on our machines!)

Second, we rename the sequences using a commmon sequence identifier:

Locality|Date|Country|Sample\_Number

In this format:
* Locality is generally the city where the sample came from. This is usually taken from the [locality](https://developers.google.com/maps/documentation/geocoding/intro#Types) entry from the Google Maps API. We chose this as it is a uniform way to ensure every sample has a locality
* Date is an eight digit number in the format YYYYMMDD. If the date is not known, we use 00000000.
* Country is the name of the country
* Sample number is to disambiguate different samples from the same city on the same date.

Note also that both locality and country match the regular expression \w+ (i.e. only contain the characters a-z, A-Z, 0-9 and \_). 

Next, we blast the sequences against the crAssphage sequence. This allows us to do two things: 1) we reverse complement any sequences that are in the wrong orientation and 2) we remove any spurious sequences that do not appear to be similar to crAssphage. We have a few of these that have leaked through our pipeline to create strains, and we are not sure where they are from.

Now that we have all the good sequences in the same direction, we align them using [muscle](www.drive5.com/muscle/muscle.html). There is a script for doing that, and if you are using the cluster, you can make it go very fast by making a virtual parallel environment and running the alignment in that. For example, to submit to the cluster and use 150 cores as an environment, you can use the command `qsub -cwd -o sge_out -e sge_err -pe make 150 ./muscle.sh`.

Once the alignment is complete we trim to remove any column that is <90% _informative characters_ and subsequently to remove any sequence that is <80%_ informative characters_. By _informative characters_, we basically mean not a hyphen (_-_). 

Finally, we use either [FastTree](microbesonline.org/fasttree/) (well, actually FastTreeDbl) or [iq-tree](http://www.iqtree.org/) to build the tree from the alignment.

As noted above, all of this is in the Makefile, however we typically run it in two steps:

```
make reverse_complement
qsub -cwd -o sge_out -e sge_err -pe make 150 ./muscle.sh
```

The first step makes the fasta file, blasts it, and reverse complements as needed.

The second step makes an alignment of the sequences using muscle.

## Making an IQtree

We also have an approach to use [iq-tree](http://www.iqtree.org/). The code is in  `iqtree.sh`. There are several issues with iq-tree. First, it checks the machine it is running on to see how many possible cores it can use. However, if you run it in a parallel environment like the command above for muscle, and force it to use more cores, it complains that you are asking it to use more cores than are available, even though you are not. Second, if you run it without specifying how many threads to use, it will determine what it thinks is the right solution, and in my case runs with two threads. However, empirical testing shows that it runs faster with more threads. If you specify the number of threads (as we do in `iqtree.sh`) at _every_ iteration it complains that it thinks it is using too many threads. Even though it runs faster. Finally, iq-tree unilaterally renames sequences, replacing | with \_. There is no reason to do this, as far as I can tell.


The files that iq-tree outputs are: (these examples are linked to PrimerC, but there are equivilent files for primers A and B).
IQ-TREE report: [seqs.C.rc.trim.aln.iqtree](PrimerC/seqs.C.rc.trim.aln.iqtree)
Maximum-likelihood tree: [seqs.C.rc.trim.aln.treefile](PrimerC/seqs.C.rc.trim.aln.treefile)
Likelihood distances: [seqs.C.rc.trim.aln.mldist](PrimerC/seqs.C.rc.trim.aln.mldist.gz)
Screen log file: [seqs.C.rc.trim.aln.log](seqs.C.rc.trim.aln.log)

Note that we store seqs.C.rc.trim.aln.mldist as a gzip compressed file in the repo to save space.


## Making the world maps

We use [cartopy](http://scitools.org.uk/cartopy/) to make the world maps. The maps are made in two stages. First we need to convert the tree to a cophenetic matrix, and then we use the matrix to find the closest elements. We need to make the matrix otherwise you need to iterate over the tree to find closest neighbors.

```
python3 ../../../bin/tree_to_cophenetic_matrix.py -t seqs.A.rc.trim.tree > seqs.A.rc.trim.matrix
python3 ../../../bin/crAssphage_cophenetic.py -i id.A.map -m seqs.A.rc.trim.matrix -o PrimerA_map.svg
```
