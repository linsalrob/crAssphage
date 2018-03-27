my $p = shift || die "primer?";


open(IN, "metadata.$p.tsv") || die "$! metadata.$p.tsv";
my @h;
while (<IN>) {
	chomp;
	my @a=split /\t/;
	unless (@h) {@h=@a; next}
	my $s = join(" ", map {"[$h[$_]=$a[$_]]"} (1 .. 9, 11, 12));
	$data{$a[10]}=">$a[0] $s";
}

open(IN, "norwegian_sanger_sequencing_crass_$p.fna") || die "$! norwegian_sanger_sequencing_crass_$p.fna";
while (<IN>) {
	if (s/^>//) {
		chomp;
		if ($data{$_}) {
			print "$data{$_}\n";
		} else {
			print STDERR "No data for $_\n";
			print ">$_\n";
		}
	} else  {
		print;
	}
}

	
