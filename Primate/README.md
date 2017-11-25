# Non-Human Primate Sequences

We identified crAssphage in non human primate sequences. Initially we had a few hits by mapping reads to crAssphage using bowtie2, but the sequences are so divergent we could not identify regions that we were interested in.

Instead, we assembled the sequences using [SPAdes](http://cab.spbu.ru/software/spades/), and then identified contigs that looked like crAssphage from amino-acid level homology. We identified ORFs in those contigs, and collected all the ORFs that are common to crAssphage and can be found in all the non-human primate genomes.

# Data:

- the Primate genome sequences (length20k.fna)
- the muscle-aligned "unambiguous core" protein sequences (present once in every primate) that were concatenated for the primate tree
- the iq-tree primate tree incl bootstraps


We used those to build a tree using iq-tree. To generate the bootstraps we used

```
iqtree-omp -s all_primates.muscle -nt 10 -alrt 1000 -bb 1000
```


### Thanks to our collaborators for sharing this data prior to publication. 
