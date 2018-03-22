
use strict;
my %count;

foreach my $f (@ARGV) {
	open(IN, $f) || die "$! $f";
	my $header = <IN>;
	my @a=split /\t/, $header;
	my ($country, $lat, $lon)=(undef, undef, undef);
	map {
		if ($a[$_] eq "country") {$country=$_}
		if ($a[$_] eq "latitude") {$lat=$_}
		if ($a[$_] eq "longitude") {$lon=$_}
		} (0 .. $#a);
	unless ($country && $lat && $lon) {die "Can't find one of | $country | $lat | $lon |"}
	while (<IN>) {
		chomp;
		my @a=split /\t/;
		$count{$a[$country]}->{$a[$lat]}->{$a[$lon]}->{$f}++;
	}

}


foreach my $c (keys %count) {
	foreach my $lat (keys %{$count{$c}}) {
		foreach my $lon (keys %{$count{$c}->{$lat}}) {
			print join("\t", $c, $lat, $lon);
			foreach my $f (@ARGV) {
				if ($count{$c}->{$lat}->{$lon}->{$f}) {
					print "\t", $count{$c}->{$lat}->{$lon}->{$f};
				}
				else {
					print "\t0";
				}
			}
			print "\n";
		}
	}
}

