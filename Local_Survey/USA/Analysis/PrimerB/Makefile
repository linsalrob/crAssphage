# We're using make to generate alignments and trees because, well why not!

PRIMER = B

help:
	@echo "Please use \`make <target>' where <target> is one of"
	@echo "\tsite_date\tcreate a tree and rename the leaves using both the site and date tags"
	@echo "\tcountry_tree\tcreate a tree and rename the leaves using the country= tag in the fasta def line"
	@echo "\tsite_tree\tcreate a tree and rename the leaves using the site= tag in the fasta def line"
	@echo "\tdate_tree\tcreate a tree and rename the leaves using the date= tag in the fasta def line"
	@echo "\ttree\t\tmake a tree named with the label of the sequences"
	@echo "\tblastn\t\tblast all the sequences against crAssphage (e.g. to check orientation)"
	@echo "\talign\t\talign the sequences using the default aligner (mafft)"
	@echo "\tseqs.fa\t\tassemble the DNA sequences into a single fasta file"
	@echo "\trenum\t\trenumber the fasta file and create a separate id.map"
	@echo "\tclustalw\trun a clustalw alignment on the sequences"
	@echo "\tmafft\t\trun a mafft alignment on the sequences (this makes a phylip format output used for building the trees)"
	@echo "\tmafft_fa\trun a mafft alignment but write the alignment in fasta format"
	@echo "\tmafft_fa_name\trun a mafft alignment but write the alignment in fasta format and rename the sequences"
	@echo ""
	@echo "You probably want to start with country_tree, as it will create the sequences, build an alignment, rename the sequences, etc etc"

site_date: seqs.$(PRIMER).tree
	python ../../../../bin/rename_trees.py -t seqs.$(PRIMER).tree -i id.$(PRIMER).map -n site -n date > seqs.$(PRIMER).sitedate.tree

date_tree: seqs.$(PRIMER).tree
	python ../../../../bin/rename_trees.py -t seqs.$(PRIMER).tree -i id.$(PRIMER).map -n date > seqs.$(PRIMER).date.tree

site_tree: seqs.$(PRIMER).tree
	python ../../../../bin/rename_trees.py -t seqs.$(PRIMER).tree -i id.$(PRIMER).map -n site > seqs.$(PRIMER).site.tree

country_tree: seqs.$(PRIMER).tree
	python ../../../../bin/rename_trees.py -t seqs.$(PRIMER).tree -i id.$(PRIMER).map -n country  > seqs.$(PRIMER).tree

tree: seqs.$(PRIMER).tree
	python ../../../../bin/rename_trees.py -t seqs.$(PRIMER).tree -i id.$(PRIMER).map > seqs.$(PRIMER).named.tree

seqs.$(PRIMER).tree: seqs.$(PRIMER).trim.aln
	FastTreeDbl -nt seqs.$(PRIMER).trim.aln > seqs.$(PRIMER).tree

blastn: seqs.$(PRIMER).fa
	blastn -db ../../../../data/JQ995537.fna -query seqs.$(PRIMER).fa -out seqs.crassphage.blastn -outfmt '6 std qlen slen'

seqs.$(PRIMER).trim.aln: seqs.$(PRIMER).aln 
	python ../../../../bin/trim_fasta_alignment.py -f seqs.$(PRIMER).aln -c 0.9 > seqs.$(PRIMER).trim.aln

align: mafft

mafft: seqs.$(PRIMER).aln

seqs.$(PRIMER).aln: seqs.$(PRIMER).fa id.$(PRIMER).map metadata.$(PRIMER).tsv
	mafft --adjustdirectionaccurately --localpair  --maxiterate 16 --reorder seqs.$(PRIMER).renum.fa > seqs.$(PRIMER).aln

clustalw: seqs.$(PRIMER).fa
	clustalw seqs.$(PRIMER).fa -OUTFILE=seqs.$(PRIMER).clustalw

metadata.$(PRIMER).tsv: id.$(PRIMER).map
	python ../../../bin/metadata.py -i id.$(PRIMER).map -o metadata.$(PRIMER).tsv
	
renum: id.$(PRIMER).map

id.$(PRIMER).map: seqs.$(PRIMER).fa
	python3 ../../../bin/rename_fasta.py -f seqs.$(PRIMER).fa -o seqs.$(PRIMER).renum.fa -i id.$(PRIMER).map

seqs.$(PRIMER).fa:
	@echo "MAKING A NEW FASTA"
	cat ../../Sequences/Primer$(PRIMER)/* > seqs.$(PRIMER).fa

clean:
	rm -f *.aln *.dnd seqs.*.fa seqs.?.renum.fa id.*.map seqs.*.dnadist  seqs.*.tree  seqs.*.tree seqs.aln.*.full seqs.*.faln  seqs.?.named.tree seqs.?.nonneg.tree
	rm -f seqs.$(PRIMER).sitedate.tree seqs.$(PRIMER).date.tree seqs.$(PRIMER).site.tree seqs.crassphage.blastn metadata.$(PRIMER).tsv
