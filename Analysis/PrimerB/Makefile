# We're using make to generate alignments and trees because, well why not!


help:
	@echo "Please use \`make <target>' where <target> is one of"
	@echo "    country_tree    create a neighbor joining tree and rename the leaves using the country= tag in the fasta def line"
	@echo "    blastn    blast all the sequences against crAssphage (e.g. to check orientation)"
	@echo "    align     align the sequences using the default aligner (mafft)"
	@echo "	   seqs.fa   assemble the DNA sequences into a single fasta file"
	@echo "    renum     renumber the fasta file and create a separate id.map"
	@echo "    clustalw  run a clustalw alignment on the sequences"
	@echo "    mafft     run a mafft alignment on the sequences"
	@echo "    neighbor_tree   create a neighbor joining tree"
	@echo ""
	@echo "You probably want to start with country_tree, as it will create the sequences, build an alignment, rename the sequences, etc etc"

country_tree: neighbor_tree
	python ../bin/rename_trees.py -t seqs.neighbor -i id.map -n country  > seqs.tree

tree: neighbor_tree
	python ../bin/rename_trees.py -t seqs.neighbor -i id.map > seqs.tree

neighbor_tree: align
	cp seqs.aln infile
	rm -f outfile 
	echo "y\n" | phylip dnadist
	mv -f outfile infile
	echo "y\n" | phylip neighbor
	mv outtree seqs.neighbor
	mv infile seqs.dnadist
	rm -f outfile


blastn: seqs.fa
	blastn -db ../data/JQ995537.fna -query seqs.fa -out seqs.crassphage.blastn -outfmt '6 std qlen slen'


align: mafft

mafft: seqs.fa renum
	mafft --adjustdirectionaccurately --localpair  --maxiterate 16 --phylipout --reorder seqs.renum.fa > seqs.aln

clustalw:
	clustalw seqs.fa

renum: seqs.fa
	perl ../bin/renumber_fasta.pl seqs.fa seqs.renum.fa id.map

seqs.fa:
	cat ../../Sequences/PrimerB/*.fasta > seqs.fa

clean:
	rm -f *.aln *.dnd seqs.fa seqs.renum.fa id.map seqs.dnadist  seqs.neighbor  seqs.tree
