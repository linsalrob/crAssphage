
use strict;
use Rob;
my $rob = new Rob;

foreach my $f (@ARGV) {
	my $fa = $rob->read_fasta($f);
	$f =~ /^(\d+)\_primer(.)\_(...)\.seq/;
	my ($seq, $prim, $dir) = ($1, $2, $3);
	unless ($seq && $prim && $dir) {
		die "Can't parse $f";
	}
	my $country = 'Israel';
	my ($lat, $lon, $city);
	if ($seq == 9382) {$lat = 31.85; $lon = 35.45; $city = "Jericho"}
	if ($seq == 9386) {$lat = 32.08; $lon = 35.18; $city = 'Salfit'}
	if ($seq == 9387) {$lat = 31.90; $lon = 35.19; $city = 'Ramallah'}
	if ($seq == 9388) {$lat = 32.18; $lon = 34.97; $city = 'Qalqilyah'}
	
	my $fat = ">Qimron_${city}_${prim}_20151202_$dir [name=Qimron lab primer $prim 20151202 $dir] [latitude=$lat] [longitude=$lon] [note=Collected in $city] [country=$country]";
	my @ids = keys %$fa;
	if ($#ids > 0) {die "More than one sequence in $f"}
	my $of = $f;
	$of =~ s/.seq/.fasta/;
	$of =~ s/^\d+/Qimron_$city/;
	if (-e $of) {die "not overwriting $of"}
	open(OUT, ">$of") || die "$! $of";
	print OUT "$fat\n", $fa->{$ids[0]}, "\n";
	close OUT;
}
