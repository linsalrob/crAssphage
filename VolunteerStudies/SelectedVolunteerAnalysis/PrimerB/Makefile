# We're using make to generate alignments and trees because, well why not!

PRIMER = B

.PHONY: clean help reverse_complement reverse_complement_tree mafft align

help:
	@echo "Please use \`make <target>' where <target> is one of"
	@echo "\tvol_date\tcreate a tree and rename the leaves using both the volunteer and date tags"
	@echo "\tdate_tree\tcreate a iq-tree and rename the leaves using the date= tag in the fasta def line"
	@echo "\ttree\t\tmake a tree named with the label of the sequences"
	@echo "\treverse_complement_tree\tReverse complement sequences as below and make an alignment and tree of them"
	@echo "\treverse_complement\tReverse complement sequences as required and filter to only those that match crassphage"
	@echo "\tblastn\t\tblast all the sequences against crAssphage (e.g. to check orientation)"
	@echo "\talign\t\talign the sequences using the default aligner (mafft)"
	@echo "\tseqs.fa\t\tassemble the DNA sequences into a single fasta file"
	@echo "\trenum\t\trenumber the fasta file and create a separate id.map"
	@echo "\tclustalw\trun a clustalw alignment on the sequences"
	@echo "\tmafft\t\trun a mafft alignment on the sequences (this makes a phylip format output used for building the trees)"
	@echo "\tmafft_fa\trun a mafft alignment but write the alignment in fasta format"
	@echo "\tmafft_fa_name\trun a mafft alignment but write the alignment in fasta format and rename the sequences"
	@echo "\treverse_complement\trenumber the sequences, run a blast to match them, and reverse complement ones in the wrong orientation"
	@echo ""
	@echo "You probably want to start with country_tree, as it will create the sequences, build an alignment, rename the sequences, etc etc"

vol_date: seqs.$(PRIMER).rc.trim.aln.treefile
	python3 ../../../bin/rename_trees.py -t seqs.$(PRIMER).rc.trim.aln.treefile -i id.$(PRIMER).map -n volunteer -n date > seqs.$(PRIMER).volunteerdate.tree

date_tree: seqs.$(PRIMER).rc.trim.aln.treefile
	python3 ../../../bin/rename_trees.py -t seqs.$(PRIMER).rc.trim.aln.treefile -i id.$(PRIMER).map -n date > seqs.$(PRIMER).date.tree

tree: seqs.$(PRIMER).rc.trim.aln.treefile
	python3 ../../../bin/rename_trees.py -t seqs.$(PRIMER).rc.trim.aln.treefile -i id.$(PRIMER).map > seqs.$(PRIMER).named.tree

seqs.$(PRIMER).rc.trim.aln.treefile: seqs.$(PRIMER).rc.trim.tree
	iq-tree -nt AUTO -s seqs.$(PRIMER).rc.trim.aln

FastTree: seqs.$(PRIMER).rc.trim.tree
	FastTreeDbl -nt seqs.$(PRIMER).rc.trim.aln > seqs.$(PRIMER).fasttree

blastn: seqs.$(PRIMER).fa
	blastn -task blastn -db ../../../data/JQ995537.fna -query seqs.$(PRIMER).fa -out seqs.crassphage.blastn -outfmt '6 std qlen slen' -evalue 0.001


seqs.$(PRIMER).trim.aln: seqs.$(PRIMER).aln
	python3 ../../../bin/trim_alignment.py -p seqs.$(PRIMER).aln -c 0.9 > seqs.$(PRIMER).trim.aln

reverse_complement_tree: seqs.$(PRIMER).rc.trim.tree

seqs.$(PRIMER).rc.trim.tree: seqs.$(PRIMER).renum.rc.fa
	muscle -in seqs.$(PRIMER).renum.rc.fa -out seqs.$(PRIMER).renum.rc.aln -maxiters 2 -diags
	python3 ../../../bin/trim_fasta_alignment.py -f seqs.$(PRIMER).renum.rc.aln -c 0.9 -r 0.8 -v > seqs.$(PRIMER).rc.trim.aln

reverse_complement: seqs.$(PRIMER).renum.rc.fa

seqs.$(PRIMER).renum.rc.fa: metadata.$(PRIMER).tsv
	blastn -task blastn -db ../../../data/JQ995537.fna -query seqs.$(PRIMER).renum.fa -out seqs.$(PRIMER).crassphage.blastn -outfmt '6 std qlen slen' -evalue 0.001
	python3 ../../../bin/reverse_complement.py -f seqs.$(PRIMER).renum.fa -b seqs.$(PRIMER).crassphage.blastn -o seqs.$(PRIMER).renum.rc.fa

metadata.$(PRIMER).tsv: id.$(PRIMER).map
	python3 ../../../bin/metadata.py -i id.$(PRIMER).map -o metadata.$(PRIMER).tsv

renum: id.$(PRIMER).map

id.$(PRIMER).map: seqs.$(PRIMER).fa
	python3 ../../../bin/rename_fasta.py -f seqs.$(PRIMER).fa -o seqs.$(PRIMER).renum.fa -i id.$(PRIMER).map

seqs.$(PRIMER).fa:
	for VOLUNTEER in D E ; do \
		echo "Reading ../../Sequences/Primer$(PRIMER)/Volunteer$$VOLUNTEER.*.fasta"; \
		cat ../../Sequences/Primer$(PRIMER)/Volunteer$$VOLUNTEER.*.fasta >> seqs.$(PRIMER).fa; \
	done


clean:
	rm -f *.aln *.dnd seqs.*.fa seqs.?.renum.fa id.*.map seqs.*.dnadist  seqs.*.tree  seqs.*.tree seqs.aln.*.full seqs.*.faln  seqs.?.named.tree seqs.?.nonneg.tree
	rm -f seqs.$(PRIMER).sitedate.tree seqs.$(PRIMER).date.tree seqs.$(PRIMER).site.tree seqs.$(PRIMER).crassphage.blastn seqs.$(PRIMER).renum.rc.fa metadata.$(PRIMER).tsv
	rm -f seqs.$(PRIMER).rc.trim.fasttree
	rm -f seqs.$(PRIMER).rc.trim.aln.bionj seqs.$(PRIMER).rc.trim.aln.iqtree seqs.$(PRIMER).rc.trim.aln.mldist seqs.$(PRIMER).rc.trim.aln.treefile seqs.$(PRIMER).rc.trim.aln.ckp.gz
	rm -f seqs.$(PRIMER).rc.trim.aln.log seqs.$(PRIMER).rc.trim.aln.model seqs.$(PRIMER).rc.trim.aln.uniqueseq.phy sgeout iqout iqerr sgeerr
