#_perl_
#

use strict;

use Rob;
my $rob = new Rob;
opendir(DIR, ".") || die "$! .";
my @files = grep {/fasta$/} readdir(DIR);
closedir(DIR);

foreach my $f (@files) {
	$f =~ m/JB.Sample(\d)([ABC])-Primer/;
	unless ($1 && $2) {print STDERR "Can't parse $f\n"; next}
	my $fa = $rob->read_fasta($f);
	open(OUT, ">>Sample$1-Primer$2.fna") || die "$! >>Sample$1-Primer$2.fna";
	foreach my $id (keys %$fa) {
		my $seq = substr($fa->{$id}, 0, 810);
		$seq = substr($seq, 10);
		#if ($f =~ /Rw/) {$seq=$rob->rc($seq)}
		print OUT ">$id\n$seq\n";
	}
	close OUT;
}

