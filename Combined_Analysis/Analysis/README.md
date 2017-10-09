# Trees and alignments of all the data.

The trees and alignments are made in several steps, although all the steps are listed in the [Makefile](PrimerA/Makefile)

First, we collect all the sequenes in fasta format. Hopefully, those sequences are symbolic links to the original data directories (they are on our machines!)

Second, we renumber the sequences as usual, so we can separate out the ID and the metadata.

Next, we blast the sequences against the crAssphage sequence. This allows us to do two things: 1) we reverse complement any sequences that are in the wrong orientation and 2) we remove any spurious sequences that do not appear to be similar to crAssphage. We have a few of these that have leaked through our pipeline to create strains, and we are not sure where they are from.

Now that we have all the good sequencces in the same direction, we align them using [muscle](www.drive5.com/muscle/muscle.html). There is a script for doing that, and if you are using the cluster, you can make it go very fast by making a virtual parallel environment and running the alignment in that. For example, to submit to the cluster and use 150 cores as an environment, you can use the command `qsub -cwd -o sge_out -e sge_err -pe make 150 ./muscle.sh`.

Once the alignment is complete we trim to remove any column that is <90% _informative characters_ and subsequently to remove any sequence that is <80%_ informative characters_. By _informative characters_, we basically mean not a hyphen (_-_). 

Finally, we use [FastTree](microbesonline.org/fasttree/) (well, actually FastTreeDbl) to build the tree from the alignment.

As noted above, all of this is in the Makefile, however we typically run it in two steps:

```
make reverse_complement
qsub -cwd -o sge_out -e sge_err -pe make 150 ./muscle.sh
```

The first step makes the fasta file, blasts it, and reverse complements as needed.

The second step does all the rest.
