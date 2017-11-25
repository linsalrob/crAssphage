# Installing the crAssphage analysis package

As you might expect, we have several dependencies that we rely on for this analysis. You will need to install these to run the phylogenetics analysis on the sequences. (However, you can run your own pipelines on the data without installing anything here - just download the sequences and start analysing them!)

## Phylogenetics

Install all of:
- [phylip](http://evolution.genetics.washington.edu/phylip/)
- [FastTreeDbl](http://microbesonline.org/fasttree/) (Note, [as described here by Aaron Darling](http://darlinglab.org/blog/2015/03/23/not-so-fast-fasttree.html), be sure to use fast tree double, not regular fast tree).
- [iq-tree](http://www.iqtree.org/)

## Similarities
- [blast+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)

## Sequence alignment

Install all three of:
- [mafft](http://mafft.cbrc.jp/alignment/software/)
- [clustalw](http://www.clustal.org/clustal2/)
- [clustalx](http://www.clustal.org/clustal2/)


## Python and Perl

We use python version 3 and perl 5 for all the other code.

### Python libraries we depend on

In addition to the usual scientific stack that you probably already have installed (numpy, matplotlib, scipy), we also use these libraries for different aspects of the analysis:

- [ete3](http://etetoolkit.org/) to parse tree files, calculate distances, and relabel leaves
- [cartopy](http://scitools.org.uk/cartopy/) to make the map figures
