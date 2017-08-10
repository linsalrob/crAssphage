
use strict;

my %data; my %seen;
while (<DATA>) {
	chomp;
	my ($file, $id1, $id, $when, $d1, $date)=split /\t/;
	my $tag = "${id}_${date}";
	$seen{$tag}++;
	$tag .= "_$seen{$tag}";
	$data{$file}=">$tag [name=EdwardsLab Volunteer Study sample $id] [date=$date] [samplecount=$seen{$tag}] [latitude=32.777166] [longitude=-117.071116= [location=San Diego, CA] [type=Fecal Samples] [country=USA]";
}

foreach my $dir ("../raw_data/") {
	opendir(DIR, $dir) || die "$! $dir";
	foreach my $f (grep {/\.seq$/} readdir(DIR)) {
		unless ($data{$f}) {print STDERR "Can't find data for $f\n"; next}
		print "$data{$f}\n";
		open(IN, "$dir/$f") || die "$! $dir/$f";
		while (<IN>) {chomp; print}
		print "\n";
	}
}




__DATA__
UB_019_D1V11F_V11F_A04.seq	D1	D	pre	4/21/17	20170421
UB_020_D1V11R_V11R_B04.seq	D1	D	pre	4/21/17	20170421
UB_017_E1V11F_V11F_G03.seq	E1	E	pre	4/21/17	20170421
UB_018_E1V11R_V11R_H03.seq	E1	E	pre	4/21/17	20170421
UB_011_MiV11F_V11F_A03.seq	Mi	Mi	pre	5/17/17	20170517
UB_012_MiV11R_V11R_B03.seq	Mi	Mi	pre	5/17/17	20170517
UB_025_Y1V11F_V11F_G04.seq	Y1	Y	pre	5/18/17	20170518
UB_026_Y1V11R_V11R_H04.seq	Y1	Y	pre	5/18/17	20170518
UB_023_WE1V11F_V11F_E04.seq	E1	E	W	5/22/17	20170522
UB_024_WE1V11R_V11R_F04.seq	E1	E	W	5/22/17	20170522
UB_013_WMiV11F_V11F_C03.seq	Mi	Mi	W	5/24/17	20170524
UB_014_WMiV11R_V11R_D03.seq	Mi	Mi	W	5/24/17	20170524
UB_005_Ke1V11F_V11F_C02.seq	Ke1	Ke	pre	5/25/17	20170525
UB_006_Ke1V11R_V11R_D02.seq	Ke1	Ke	pre	5/25/17	20170525
UB_021_WD1V11F_V11F_C04.seq	D1	D	W	6/1/17	20170601
UB_022_WD1V11R_V11R_D04.seq	D1	D	W	6/1/17	20170601
UB_009_WD2V11F_V11F_G02.seq	D2	D	W	6/1/17	20170601
UB_010_WD2V11R_V11R_H02.seq	D2	D	W	6/1/17	20170601
UB_007_WE2V11F_V11F_E02.seq	E2	E	W	6/1/17	20170601
UB_008_WE2V11R_V11R_F02.seq	E2	E	W	6/1/17	20170601
UB_015_DD3V11F_V11F_E03.seq	D3	D	D	6/8/17	20170608
UB_016_DD3V11R_V11R_F03.seq	D3	D	D	6/8/17	20170608
UB_003_WMi4V11F_V11F_A02.seq	Mi	Mi	W	6/22/17	20170622
UB_004_WMi4V11R_V11R_B02.seq	Mi	Mi	W	6/22/17	20170622
UB_001_DD5V11F_V11F_G01.seq	D5	D	D	6/23/17	20170623
UB_002_DD5V11R_V11R_H01.seq	D5	D	D	6/23/17	20170623
