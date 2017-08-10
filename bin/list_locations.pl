#__perl__
#


=pod

Generate a list of the locations of the different samples. This is used 
to help write the paper!

=cut

use strict;
use Rob;
my $rob = new Rob;
my %count;

my $dir = shift || die "Directory?";
opendir(DIR, $dir) || die $!;
foreach my $f (grep {$_ !~ /^\./} readdir(DIR)) {
	my $fa = $rob->read_fasta("$dir/$f");
	foreach my $t (keys %$fa) {
		my ($lat, $lon, $co)=("", "", "");
		if ($t =~ m/latitude=(.*?)\]/) {$lat=$1}
		else {print STDERR "$f :: NO LAT: $t\n"}
		if ($t =~ m/longitude=(.*?)\]/) {$lon=$1}
		else {print STDERR "$f :: NO LON: $t\n"}
		if ($t =~ m/country=(.*?)\]/) {$co=$1}
		else {print STDERR "$f :: NO COUNTRY: $t\n"}
		my $loc = join("_", $co, $lat,$lon);
		$count{$loc}++;
	}
}

foreach my $k (keys %count) {print "$k\t$count{$k}\n"}
		
