# We're using make to generate alignments and trees because, well why not!


help:
	@echo "Please use \`make <target>' where <target> is one of"
	@echo "    country_tree    create a neighbor joining tree and rename the leaves using the country= tag in the fasta def line"
	@echo "    blastn    blast all the sequences against crAssphage (e.g. to check orientation)"
	@echo "    align     align the sequences using the default aligner (mafft)"
	@echo "	   seqs.fa   assemble the DNA sequences into a single fasta file"
	@echo "    renum     renumber the fasta file and create a separate id.map"
	@echo "    clustalw  run a clustalw alignment on the sequences"
	@echo "    mafft     run a mafft alignment on the sequences (this makes a phylip format output used for building the trees)"
	@echo "    mafft_fa  run a mafft alignment but write the alignment in fasta format"
	@echo "    neighbor_tree   create a neighbor joining tree"
	@echo ""
	@echo "You probably want to start with country_tree, as it will create the sequences, build an alignment, rename the sequences, etc etc"

country_tree: seqs.country.tree

seqs.country.tree: seqs.neighbor
	python3 ../../../bin/rename_trees_fast.py -t seqs.neighbor -i id.map -n country  > seqs.country.tree

tree: seqs.named.tree

seqs.named.tree: seqs.neighbor
	python3 ../../../bin/rename_trees_fast.py -t seqs.neighbor -i id.map > seqs.named.tree

neighbor_tree: seqs.neighbor

seqs.neighbor: seqs.trim.aln
	FastTreeDbl -nt seqs.trim.aln > seqs.neighbor


blastn: seqs.fa
	blastn -db ../data/JQ995537.fna -query seqs.fa -out seqs.crassphage.blastn -outfmt '6 std qlen slen'


seqs.trim.aln: seqs.aln
	python3 ../../../bin/trim_alignment.py -p seqs.aln -c 0.9 > seqs.trim.aln

align: mafft

mafft_fa: seqs.faln

seqs.faln: seqs.fa renum
	mafft --adjustdirectionaccurately --localpair  --maxiterate 16 --reorder seqs.renum.fa > seqs.faln

mafft: seqs.aln

seqs.aln: seqs.fa id.map
	mafft --adjustdirectionaccurately --localpair  --maxiterate 16 --phylipout --reorder seqs.renum.fa > seqs.aln

clustalw: seqs.clustalw

seqs.clustalw:
	# clustalw -INFILE=seqs.fa -OUTFILE=seqs.clustalw
	mafft --adjustdirectionaccurately --localpair  --maxiterate 16 --reorder  --clustalout --thread 8 seqs.renum.fa > seqs.clustalw

renum: id.map

id.map: seqs.fa
	perl ../../../bin/renumber_fasta.pl seqs.fa seqs.renum.fa id.map

seqs.fa:
	#ln -s gretel_haplotypes_pcrA.fasta seqs.fa
	ln -s haplotypes_country.fasta seqs.fa

clean:
	rm -f *.aln *.dnd seqs.fa seqs.renum.fa id.map seqs.dnadist  seqs.neighbor  seqs.tree 
	rm -f seqs.aln.full seqs.faln  seqs.named.tree seqs.nonneg.tree seqs.clustalw seqs.country.tree

