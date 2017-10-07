# Data analysis

We have included our data analysis in the repository, but we also include all the raw data, and so you are completely free to perform your own data analysis if you want!

# Makefile

We're using the [GNU Make](https://www.gnu.org/software/make/) program to automatically run the analysis whenever we add new sequences. This is originally designed for compiling software, but it allows us to create some dependencies and run things sequentially. For example, you can't create a tree until you have an alignment, and you can't create an alignment until you have some sequences! 

If you look inside the [Makefile](PrimerA/Makefile), you will likely notice a lot of common commands. The advantage of using the Makefile is that we can use the same parameters and approach each time we perform an analysis.

Our general pipeline is:

1. Concatenate all the sequences.
2. Renumber all the sequences. Some of the software we use doesn't handle long sequence names, so we just number them starting at 1.
3. Use [mafft](http://mafft.cbrc.jp/alignment/software/) to align the sequences. MAFFT has several nice attributes, including the ability to automatically reverse-complement a sequence, allow local alignments rather than global alignments, and so on. 
4. We trim the sequences to informative sites, typically those that in >90% of sequences are not *-* or *N*.
5. Use dnadist and neighbor in the [phylip](http://evolution.genetics.washington.edu/phylip/) package to score the alignments and generate the trees.
6. Rename the trees using some criteria you decide.

If you have [GNU Make](https://www.gnu.org/software/make/) installed, you should be able to type:

```
make country_tree
```

and all these steps will be completed for you. The tree in `seqs.tree` will be labeled with the country of origin of the sequences. 

We suggest that you use the excellent [figtree](http://tree.bio.ed.ac.uk/software/figtree/) software to visualize the tree that results.

Of course, you are welcome to develop (and contribute) your own sequence analysis approach.
