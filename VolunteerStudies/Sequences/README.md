# Cleaning the sequences

The ab1 files were manually inspected using FinchTV from Geospiza and the base calls written in fasta format.

Forward and reverse reads for each sample were merged into a single file and aligned using mafft:

```
mafft --adjustdirectionaccurately  --maxiterate 16   --thread 8  seqs.fna > seqs.aln
```

Alignments were visually inspected using clustalx to make sure that there was enough overlap between the two reads to ensure accurate alignment.

```
clustalx seqs.aln
```

The consensus base calls from the alignment were called and written out:

```
fasta_consensus.pl seqs.aln  > ../final_sequences/seqs.fna
```

Finally, the sequences were renamed to make sure we can keep track of everything.
