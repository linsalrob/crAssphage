# We're using make to generate alignments and trees because, well why not!

PRIMER = A
VOLUNTEER = E

.PHONY: clean help reverse_complement reverse_complement_tree mafft align

help:
	@echo "Please use \`make <target>' where <target> is one of"
	@echo "\tvol_date\tcreate a tree and rename the leaves using both the volunteer and date tags"
	@echo "\tvolunteer_tree\tcreate a tree and rename the leaves using the volunteer= tag in the fasta def line"
	@echo "\tdate_tree\tcreate a tree and rename the leaves using the date= tag in the fasta def line"
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
	@echo "You probably want to start with vol_date, as it will create the sequences, build an alignment, rename the sequences, etc etc"

vol_date: seqs.$(VOLUNTEER).$(PRIMER).tree
	python3 ../../../../bin/rename_trees.py -t seqs.$(VOLUNTEER).$(PRIMER).tree -i id.$(VOLUNTEER).$(PRIMER).map -n volunteer -n date > seqs.$(VOLUNTEER).$(PRIMER).volunteerdate.tree

date_tree: seqs.$(VOLUNTEER).$(PRIMER).tree
	python3 ../../../../bin/rename_trees.py -t seqs.$(VOLUNTEER).$(PRIMER).tree -i id.$(VOLUNTEER).$(PRIMER).map -n date > seqs.$(VOLUNTEER).$(PRIMER).date.tree

vol_tree: seqs.$(VOLUNTEER).$(PRIMER).tree
	python3 ../../../../bin/rename_trees.py -t seqs.$(VOLUNTEER).$(PRIMER).tree -i id.$(VOLUNTEER).$(PRIMER).map -n volunteer > seqs.$(VOLUNTEER).$(PRIMER).volunteer.tree

tree: seqs.$(VOLUNTEER).$(PRIMER).tree
	python3 ../../../../bin/rename_trees.py -t seqs.$(VOLUNTEER).$(PRIMER).tree -i id.$(VOLUNTEER).$(PRIMER).map > seqs.$(VOLUNTEER).$(PRIMER).named.tree

seqs.$(VOLUNTEER).$(PRIMER).tree: seqs.$(VOLUNTEER).$(PRIMER).trim.aln
	FastTreeDbl -nt seqs.$(VOLUNTEER).$(PRIMER).trim.aln > seqs.$(VOLUNTEER).$(PRIMER).tree

blastn: seqs.$(VOLUNTEER).$(PRIMER).fa
	blastn -task blastn -db ../../../../data/JQ995537.fna -query seqs.$(VOLUNTEER).$(PRIMER).fa -out seqs.$(VOLUNTEER).$(PRIMER).crassphage.blastn -outfmt '6 std qlen slen' -evalue 1e-12


seqs.$(VOLUNTEER).$(PRIMER).trim.aln: seqs.$(VOLUNTEER).$(PRIMER).aln
	python3 ../../../../bin/trim_fasta_alignment.py -f seqs.$(VOLUNTEER).$(PRIMER).aln -c 0.9 > seqs.$(VOLUNTEER).$(PRIMER).trim.aln

seqs.all.aln:
	muscle -in seqs.all.renum.fa -out seqs.all.aln -maxiters 2 -diags

align: mafft

mafft: seqs.$(VOLUNTEER).$(PRIMER).aln

seqs.$(VOLUNTEER).$(PRIMER).aln: seqs.$(VOLUNTEER).$(PRIMER).fa id.$(VOLUNTEER).$(PRIMER).map metadata.$(VOLUNTEER).$(PRIMER).tsv
	mafft --adjustdirectionaccurately --localpair  --maxiterate 16 --reorder seqs.$(VOLUNTEER).$(PRIMER).renum.fa > seqs.$(VOLUNTEER).$(PRIMER).aln

phylipoutput:
	mafft --adjustdirectionaccurately --localpair  --maxiterate 16 --phylipout --reorder seqs.$(VOLUNTEER).$(PRIMER).renum.fa > seqs.$(VOLUNTEER).$(PRIMER).aln

mafft_fa: seqs.$(VOLUNTEER).$(PRIMER).fa id.$(VOLUNTEER).$(PRIMER).map
	mafft --adjustdirectionaccurately --localpair  --maxiterate 16 --reorder seqs.$(VOLUNTEER).$(PRIMER).renum.fa > seqs.$(VOLUNTEER).$(PRIMER).faln

clustalw: seqs.$(VOLUNTEER).$(PRIMER).fa
	clustalw seqs.$(VOLUNTEER).$(PRIMER).fa -OUTFILE=seqs.$(VOLUNTEER).$(PRIMER).clustalw

reverse_complement_tree: seqs.$PRIMER.rc.trim.tree

seqs.$PRIMER.rc.trim.tree: seqs.$(VOLUNTEER).$(PRIMER).renum.rc.fa
	@echo "Running muscle may take a while. You might want to run these steps on the cluster"
	muscle -in seqs.$(VOLUNTEER).$(PRIMER).renum.rc.fa -out seqs.$(VOLUNTEER).$(PRIMER).renum.rc.aln -maxiters 2 -diags
	python3 ../../../../bin/trim_fasta_alignment.py -f seqs.$(VOLUNTEER).$(PRIMER).renum.rc.aln -c 0.9 -r 0.8 -v > seqs.$(VOLUNTEER).$(PRIMER).rc.trim.aln
	fasttree -nt seqs.$(VOLUNTEER).$(PRIMER).rc.trim.aln > seqs.$(VOLUNTEER).$(PRIMER).rc.trim.tree

reverse_complement: seqs.$(VOLUNTEER).$(PRIMER).renum.rc.fa

seqs.$(VOLUNTEER).$(PRIMER).renum.rc.fa: id.$(VOLUNTEER).$(PRIMER).map metadata.$(VOLUNTEER).$(PRIMER).tsv
	blastn -task blastn -db ../../../../data/JQ995537.fna -query seqs.$(VOLUNTEER).$(PRIMER).renum.fa -out seqs.$(VOLUNTEER).$(PRIMER).crassphage.blastn -outfmt '6 std qlen slen' -evalue 0.001
	python3 ../../../../bin/reverse_complement.py -f seqs.$(VOLUNTEER).$(PRIMER).renum.fa -b seqs.$(VOLUNTEER).$(PRIMER).crassphage.blastn -o seqs.$(VOLUNTEER).$(PRIMER).renum.rc.fa

metadata.$(VOLUNTEER).$(PRIMER).tsv: id.$(VOLUNTEER).$(PRIMER).map
	python3 ../../../../bin/metadata.py -i id.$(VOLUNTEER).$(PRIMER).map -o metadata.$(VOLUNTEER).$(PRIMER).tsv

renum: id.$(VOLUNTEER).$(PRIMER).map

id.$(VOLUNTEER).$(PRIMER).map: seqs.$(VOLUNTEER).$(PRIMER).fa
	python3 ../../../../bin/rename_fasta.py -f seqs.$(VOLUNTEER).$(PRIMER).fa -o seqs.$(VOLUNTEER).$(PRIMER).renum.fa -i id.$(VOLUNTEER).$(PRIMER).map -t volunteer -t date

seqs.$(VOLUNTEER).$(PRIMER).fa:
	cat ../../../Sequences/Primer$(PRIMER)/Volunteer$(VOLUNTEER).*.fasta > seqs.$(VOLUNTEER).$(PRIMER).fa

clean:
	rm -f *.aln *.dnd seqs.*.fa seqs.?.renum.fa id.*.map seqs.*.dnadist  seqs.*.neighbor  seqs.*.tree seqs.aln.*.full seqs.*.faln  seqs.?.named.tree seqs.?.nonneg.tree
	rm -f seqs.$(VOLUNTEER).$(PRIMER).sitedate.tree seqs.$(VOLUNTEER).$(PRIMER).date.tree seqs.$(VOLUNTEER).$(PRIMER).site.tree seqs.$(VOLUNTEER).$(PRIMER).crassphage.blastn seqs.$(VOLUNTEER).$(PRIMER).renum.rc.fa metadata.$(VOLUNTEER).$(PRIMER).tsv
