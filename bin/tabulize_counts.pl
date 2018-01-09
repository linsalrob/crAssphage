#!/usr/bin/env perl
use strict;
my $count;
my $p; my $l;
while (<>) {
	chomp;
	if (m/primer=(\w+)/) {$p=$1} else {print STDERR "No primer found\n"}
	if (m/locality=(.*?)\]/) {$l=$1} else {print STDERR "No locality found\n"}
	$count->{$l}->{$p}++;
}

my @primers = ("A", "B", "C");
foreach my $l (sort {$a cmp $b} keys %$count) {
	print $l;
	foreach my $p (@primers) {
		print ",";
		if ($count->{$l}->{$p}) {
			print $count->{$l}->{$p};
		} else {
			print "0";
		}
	}
	print "\n";
}


