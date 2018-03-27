use strict;
open(IN, "metadata.tsv") || die "$! metadata.tsv";
my @h; my %data; my %date;
while (<IN>) {
	chomp;
	my @a=split /\t/;
	unless (@h) {@h=@a; next}
	$data{$a[0]} = join(" ", map {"[$h[$_]=$a[$_]]"} (0 .. $#a));
	$date{$a[0]}=$a[1];
}
close IN;


opendir(DIR, "GRETEL") || die "can't open GRETEL";

my $count=0;
foreach my $sample (grep {$_ !~ /^\./} readdir(DIR)) {
	unless ($data{$sample}) {print STDERR "No metadata for $sample\n"; next}
	foreach my $primer (qw[A B C]) {
		open(PCR, ">>pcr$primer.fasta") || die "can't append to pcr$primer.fasta";
		if (-e "GRETEL/$sample/gretel/pcr$primer.fasta") {
			open(FA, "GRETEL/$sample/gretel/pcr$primer.fasta") || die "$! GRETEL/$sample/gretel/pcr$primer.fasta";
			while (<FA>) {
				if (/^>/) {
					$count++;
					print PCR ">Oslo_${date{$sample}}_${sample}_Norway_$count [name=Cinek.Norway.metagenome.$sample.$count] [primer=$primer] $data{$sample}\n";
				} else {
					print PCR;
				}
			}
		}
		close PCR;
	}
}


