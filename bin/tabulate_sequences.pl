#__perl__
#

use strict;
use Rob;
my $rob = new Rob;

my %tags;
my %seqs;

my %primer = ("A" => "PrimerA", "B"=>"PrimerB", "C"=>"PrimerC");


foreach my $p (keys %primer) {
	opendir(DIR, $primer{$p}) || die "Can't open dir $primer{$p}";
	foreach my $f (grep {/fasta/} readdir(DIR)) {
			my $fa = $rob->read_fasta("$primer{$p}/$f");
			foreach my $id (keys %$fa) {
				my $newid = $id . " [primer=$p]";
				$seqs{$newid}=$fa->{$id};
				while ($newid =~ s/\[(\w+)=//) {$tags{$1}=1}
				if ($newid =~ /\[/) {print STDERR "Did we get all the tags in $id ::: $newid remaining\n"}
			}
		}
	}

#my @tags = sort {$a cmp $b} keys %tags;


my @tags = (qw[primer name contact university latitude latlon locality country sample_date date yearmonth site address location altitude longitude note sample sample_number sex source type]);


print STDERR join("\n", @tags, "\n");

print join("\t", "Label", @tags, "Sequence"), "\n";
foreach my $i (keys %seqs) {
	my $label = $i;
	$label =~ s/\s.*$//;
	print $label;
	foreach my $tag (@tags) {
		if ($i =~ /\[$tag=(.*?)\]/) {print "\t$1"} else {print "\t"}
	}
	print "\t" . uc($seqs{$i}) . "\n";
}
