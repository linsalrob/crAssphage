# Haplotype data

We created the strain data using [Gretel](https://github.com/linsalrob/crAssphage/tree/master/Metagenomes/Haplotypes) and from that created two files:

* [gretel_strains_pcrA.fasta](gretel_strains_pcrA.fasta) has all the strain sequences that we found. There are 15,221 sequences in this file, and so we haven't been able to build a tree for the sequences!
* [strains_country.fasta](strains_country.fasta) has those sequences filtered to the set for which we know the country of origin and the sampling date. There are only 2,853 sequences in this file so it is easier to build a tree from these.

There is an [unlabeled tree](strains_country.fasta) and associated [metadata](metadata.tsv) so you can make your own analyses. There is also a [tree labeled with the country of isolation](seqs.country.tree)
