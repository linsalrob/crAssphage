#__perl__
#

=pod

Convert a clustalw alignment to a fasta file.

=cut

use strict;

my $file = shift || die "File to convert from clustalw to fasta?";
my %seqs;
open(IN, $file) || die "$! $file";
my $l=<IN>;
unless ($l =~ /^CLUSTAL/) {die "The first line of $file does not begin CLUSTAL. Is it a clustalw file?"}
while (<IN>) {
	next if (/^\s+$/);
	next unless (/^\S/);
	chomp;
	m/^(.*?)\s+(\S+)$/;
	$seqs{$1}.=$2;
}
close IN;

map {print ">$_\n$seqs{$_}\n"} keys %seqs;

