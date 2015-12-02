# We're using make to generate alignments and trees because, well why not!



rename_tree: tree
	python /home/redwards/Dropbox/GitHubs/EdwardsLab/trees/rename_trees.py -t seqs.neighbor -i id.map > seqs.tree

tree: align
	cp seqs.aln infile
	rm -f outfile 
	echo "y\n" | phylip dnadist
	mv -f outfile infile
	echo "y\n" | phylip neighbor
	mv outtree seqs.neighbor
	mv infile seqs.dnadist
	rm -f outfile


blastn: seqs.fa
	blastn -db ~/Dropbox/Phage/crAssphage/JQ995537.fna -query seqs.fa -out seqs.crassphage.blastn -outfmt '6 std qlen slen'


align: mafft

mafft: seqs.fa renum
	mafft --localpair  --maxiterate 16 --phylipout --reorder seqs.renum.fa > seqs.aln

clustalw:
	clustalw seqs.fa

renum: seqs.fa
	perl /home/redwards/bioinformatics/bin/renumber_fasta.pl seqs.fa seqs.renum.fa id.map

seqs.fa:
	cat ../../Sequences/PrimerC/*.fasta > seqs.fa

clean:
	rm -f *.aln *.dnd seqs.fa seqs.renum.fa id.map
