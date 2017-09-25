#__perl__
#

use strict;
use Rob;
my $rob = new Rob;
# generate a consensus froma fasta file. This assumes that all sequences are the same length, eg. if they are the product of a gapped alignment
#

my $f = shift || die "fasta file?";
my $fa = $rob->read_fasta($f);
my @ids = keys %$fa;
print ">$f\n";
for (my $i=0; $i<=length($fa->{$ids[0]}); $i++) {
	my %s;
	map {
		my $base = substr($fa->{$_}, $i, 1);
		if ($base ne "-") {$s{$base}++};
	} @ids;
	my @bases = sort {$s{$b} <=> $s{$a}} keys %s;
	print "$bases[0]";
}
print "\n";
