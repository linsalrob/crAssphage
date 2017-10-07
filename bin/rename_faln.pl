#__perl__
#

use strict;

# Rename a fasta alignment based on a the id.map file
#


my $usage=<<EOF;
$0 <id.map> <seqs.faln>

id.map is the standard tsv with id number and metadata
seqs.faln is the fasta alignment file

EOF

my $idmapf = shift || die $usage;
my $falnf = shift || die $usage;

open(IN, $idmapf) || die "$! $idmapf";
my %id;
while (<IN>) {
	chomp;
	my @a=split /\t/;
	$a[1] =~ s/\s+.*$//;
	$id{$a[0]}=$a[1];
}
close IN;

open(IN, $falnf) || die "$! $falnf";
while (<IN>) {
	if (/^>/){
		m/(\d+)/;
		if ($id{$1}) {
			print ">$id{$1}\n";
		} else {
			print STDERR "NO ID FOR $1\n";
			print;
		}
	} else {
		print;
	}
}
