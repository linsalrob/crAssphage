#!/usr/bin/perl 
#
# renumber fasta sequences and write out an id map

use strict;
my $c=0;
my $f=shift || die "fasta file?";
my $o=shift || die "fasta output file";
my $i=shift || die "id map file to write";

open(IN, "$f") || die $!;
open(OUT, ">$o") || die $!;
open(ID, ">$i") || die $!;
while(<IN>) {
	if (s/^>//) {
		print ID "$c\t$_";
		print OUT ">", $c++, "\n";
	}
	else {
		print OUT;
	}
}
close IN;
close OUT;
close ID;
