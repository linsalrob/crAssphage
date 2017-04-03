# We're using make to generate alignments and trees because, well why not!



PRIMER = A


help:
	@echo "Please use \`make <target>' where <target> is one of"
	@echo "    tree            create a neighbor joining tree and rename the leaves using the id in the fasta def line"
	@echo "    country_tree    create a neighbor joining tree and rename the leaves using the country= tag in the fasta def line"
	@echo "    blastn    blast all the sequences against crAssphage (e.g. to check orientation)"
	@echo "    align     align the sequences using the default aligner (mafft)"
	@echo "    seqs.$(PRIMER).fa assemble the DNA sequences into a single fasta file"
	@echo "    renum     renumber the fasta file and create a separate id.$(PRIMER).map"
	@echo "    clustalw  run a clustalw alignment on the sequences"
	@echo "    mafft     run a mafft alignment on the sequences (this makes a phylip format output used for building the trees)"
	@echo "    mafft_fa  run a mafft alignment but write the alignment in fasta format"
	@echo "    neighbor_tree   create a neighbor joining tree"
	@echo ""
	@echo "You probably want to start with tree, as it will create the sequences, build an alignment, rename the sequences, etc etc"
	@echo ""
	@echo "You should define the primer you want to analyse with PRIMER=A, PRIMER=B, or PRIMER=C on the command line. The default is PRIMER=A"

site_tree: neighbor_tree
	python ../../Analysis/bin/rename_trees.py -t seqs.nonneg.$(PRIMER).tree -i id.$(PRIMER).map -n site > seqs.$(PRIMER).tree

country_tree: neighbor_tree
	python ../../Analysis/bin/rename_trees.py -t seqs.nonneg.$(PRIMER).tree -i id.$(PRIMER).map -n country  > seqs.$(PRIMER).tree

tree: neighbor_tree
	python ../../Analysis/bin/rename_trees.py -t seqs.nonneg.$(PRIMER).tree -i id.$(PRIMER).map > seqs.named.$(PRIMER).tree

neighbor_tree: seqs.nonneg.$(PRIMER).tree

seqs.nonneg.$(PRIMER).tree: seqs.$(PRIMER).aln
	cp seqs.$(PRIMER).aln infile
	rm -f outfile 
	echo "y\n" | phylip dnadist
	mv -f outfile infile
	echo "y\n" | phylip neighbor
	mv outtree seqs.$(PRIMER).neighbor
	mv infile seqs.$(PRIMER).dnadist
	rm -f outfile
	python ../../Analysis/bin/negative_branch_lengths.py -t seqs.$(PRIMER).neighbor > seqs.nonneg.$(PRIMER).tree


blastn: seqs.$(PRIMER).fa
	blastn -db ../../data/JQ995537.fna -query seqs.$(PRIMER).fa -out seqs.crassphage.blastn -outfmt '6 std qlen slen'


trim_seqs: seqs.$(PRIMER).aln
	
	
seqs.$(PRIMER).aln: seqs.aln
	mv seqs.$(PRIMER).aln seqs.aln.$(PRIMER).full
	python ../../Analysis/bin/trim_alignment.py -p seqs.aln.$(PRIMER).full -c 0.9 > seqs.$(PRIMER).aln

seqs.aln: align

align: mafft

mafft_fa: seqs.$(PRIMER).fa renum
	mafft --adjustdirectionaccurately --localpair  --maxiterate 16 --reorder seqs.renum.$(PRIMER).fa > seqs.$(PRIMER).faln

mafft: seqs.$(PRIMER).fa renum
	mafft --adjustdirectionaccurately --localpair  --maxiterate 16 --phylipout --reorder seqs.renum.$(PRIMER).fa > seqs.$(PRIMER).aln

clustalw: seqs.$(PRIMER).fa
	clustalw seqs.$(PRIMER).fa -OUTFILE=seqs.$(PRIMER).clustalw

renum: id.$(PRIMER).map

id.$(PRIMER).map: seqs.$(PRIMER).fa
	perl ../../Analysis/bin/renumber_fasta.pl seqs.$(PRIMER).fa seqs.renum.$(PRIMER).fa id.$(PRIMER).map

seqs.$(PRIMER).fa:
	@echo "MAKING A NEW FASTA"
	cat Edwards_$(PRIMER)*.fasta > seqs.$(PRIMER).fa

clean:
	rm -f *.aln *.dnd seqs.*.fa seqs.renum.*.fa id.*.map seqs.*.dnadist  seqs.*.neighbor  seqs.*.tree seqs.aln.*.full seqs.*.faln  seqs.named.*.tree seqs.nonneg.*.tree
