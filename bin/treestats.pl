#!/usr/bin/perl -w
use strict;
use Data::Dumper;

=pod

Usage: perl treestats.pl Analysis/PrimerA/seqs.neighbor Analysis/PrimerA/id.map 10000 name longitude country

The output for me is:

perl treestats.pl Analysis/PrimerA/seqs.neighbor Analysis/PrimerA/id.map 10000 name longitude country
 # metadata       stat    pvalue  permutations
 # country 17      0       10000
 # longitude       8       0.0001  10000
 # name    2       0.0048  10000

I needed to change two additional things in the metadata file Analysis/PrimerA/id.map:
- The metadata fields other than "name", "longitude", and "country" need to be completed in all lines

What the script does is it takes the tree, and scores how well the indicated phenotype (e.g. "country") clusters on the branches. The statistic is the number of 
branching partitions where all leaves have an identical value for the given metadata field - this value is in the second column of the output (e.g. 17 branches
exist in the tree where all leaves are from the same country).

Then it shuffles the labels of the tree N times (in this case N=10000), each time re-calculates the clustering score, and outputs the fraction of times that
the score in the randomized tree is equal or higher than the score in the real tree (e.g. 48 times for the "name" field). This fraction is equal to the P-value
that the clustering is significant - so all three metadata fields "name", "longitude", and "country" are significantly clustered in the tree.

The good thing about this approach is that it retains the shape of the tree (branch order etc) which may be important for the statistics.

=cut


my $usage = "Usage:\t$0 treefile.ph metadata_file Pvalue_permutations metadata1 metadata2 ...\n";
my $nrand = 0;
if (scalar @ARGV <= 3) {
	die $usage; }
for (my $i = 0; $i < 2; ++$i) {
	if (! (-e $ARGV[$i])) {
		die "No such file $ARGV[$i]: $!\n$usage"; } }
$nrand = $ARGV[2];
my %metadata_todo = ();
for (my $i = 3; $i < scalar @ARGV; ++$i) {
	++$metadata_todo{$ARGV[$i]}; }

sub shuffle (\@);

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#	READ BRACKETNOTATION FROM TREEFILE
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
my $bracketnotation = "";
open (IN, "<$ARGV[0]") or die "Can't open file $ARGV[0]: $!\n";
while (my $line = <IN>) {
	chomp ($line);
	$line =~ s/_R_//g; # this removes the notation from mafft that the reverse complement was used
	$bracketnotation .= $line; }
close (IN);
$bracketnotation =~ s/\s+/ /g;
my $bracketnotation_size = length ($bracketnotation);

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#	DETERMINE ALL PARTITIONS IN BRACKETNOTATION
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
my %partition = ();
for (my $position_rightbracket = 0; $position_rightbracket < $bracketnotation_size; ++$position_rightbracket) {
	if (substr ($bracketnotation, $position_rightbracket, 1) eq ")") {
		BACKTRACE: for (my $position_leftbracket = $position_rightbracket - 1; $position_leftbracket >= 0; --$position_leftbracket) {
			if (substr ($bracketnotation, $position_leftbracket, 1) eq "(") {
				my $leftbracket_already_used = 0;
				foreach my $used_position_rightbracket (sort keys %partition) {
					if ($partition{$used_position_rightbracket} == $position_leftbracket) {
						++$leftbracket_already_used; } }
				if ($leftbracket_already_used == 0) {
					$partition{$position_rightbracket} = $position_leftbracket;
					last BACKTRACE; } } } } }
foreach my $position_rightbracket (sort keys %partition) {
	$partition{$position_rightbracket} = substr ($bracketnotation, $partition{$position_rightbracket}, $position_rightbracket - $partition{$position_rightbracket} + 1); }


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#	READ METADATA FROM FILE
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
my %metadata = ();
my %cols = ();
open (IN, "<$ARGV[1]") or die "Can't open file $ARGV[1]: $!\n";
while (my $line = <IN>) {
	chomp ($line);
	my @this = split /\t/, $line;

# UNCOMMENT THIS PART IF YOU WANT TO USE A TAB-DELIMITED TEXT FILE FOR METADATA
#	if (scalar keys %cols == 0) {
#		foreach my $v (@this) {
#			++$cols{$v}; } }
#	else {
#		for (my $i = 1; $i < scalar @cols; ++$i) {
#			$metadata{$this[0]}{$cols[$i]} = $this[$i]; } } }

# UNCOMMENT THIS PART IF YOU WANT TO USE THE OTHER FILE FOR METADATA
	my @fields = split / \[/, $this[1];
	for (my $i = 1; $i < scalar @fields; ++$i) {
		$fields[$i] =~ s/\]\s*\Z//;
		my ($k, $v) = split /=/, $fields[$i];
		++$cols{$k};
		$metadata{$this[0]}{$k} = $v; }
	foreach my $md (keys %metadata_todo) {
		if (! (exists ($metadata{$this[0]}{$md}))) {
			die "Error:\tmetadata field '$md' not available for node $this[0] in $ARGV[1]\n$usage"; } } }

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#	DETERMINE SPECIES IN ALL PARTITIONS
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
my %species_in_clades = ();
foreach my $rb (sort { $a <=> $b; } keys %partition) {
	my $tmp_species = "$partition{$rb}";
	$tmp_species =~ s/[\):][\-\w\.]+//g;
	my @species_file = split /[\(\),;]+/, $tmp_species;
	@species_file = sort (@species_file);
	while ($species_file[0] eq "") {
		shift (@species_file); }
	for (my $j = 0; $j < scalar @species_file; ++$j) {
		++$species_in_clades{$rb}{$species_file[$j]}; } }

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#	LIST ALL SPECIES
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
my @partitions_by_size = sort { length ($partition{$a}) <=> length ($partition{$b}); } keys %partition;
my %species = ();
my @species_list = ();
foreach my $sp (keys %{$species_in_clades{$partitions_by_size[-1]}}) {
	$species{$sp} = $sp;
	push (@species_list, $sp); }

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#	COUNT METADATA VALUES FOR ALL CLADES
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
my %scores = ();
# calculate scores nrand times; for the first (shuf=0) the species are not randomized and for the rest they are
for (my $shuf = 0; $shuf <= $nrand; ++$shuf) {
	# which metadata fields to test?
	foreach my $md (keys %metadata_todo) {
		$scores{$md}[$shuf] = 0;
		# calculate a value for all partitions rb
		foreach my $rb (keys %partition) {
			my %rb_values = ();
			my $n = 0;
			foreach my $sp (keys %{$species_in_clades{$rb}}) {
				if (!defined $metadata{$species{$sp}}{$md}) {
					die "No metadata $md defined for species $sp";
				}
				++$rb_values{$metadata{$species{$sp}}{$md}};
				++$n; }
			# count if this is a perfect partition
			my @rb_values_sorted = sort { $rb_values{$a} <=> $rb_values{$b}; } keys %rb_values;
			if ($rb_values{$rb_values_sorted[-1]} == $n) {
				++$scores{$md}[$shuf]; } } }

	# after first calculating the values for the unshuffled tree, shuffle the tree
	@species_list = shuffle (@species_list);

	my $pos = 0;
	foreach my $sp (keys %species) {
		$species{$sp} = $species_list[$pos];
		++$pos; } }

print "#metadata\tstat\tpvalue\tpermutations\n";
foreach my $md (keys %scores) {
	my $pvalue = 0;
	for (my $shuf = 1; $shuf < scalar @{$scores{$md}}; ++$shuf) {
		if ($scores{$md}[$shuf] >= $scores{$md}[0]) {
			++$pvalue; } }
	$pvalue = $pvalue / $nrand;
	print "$md\t$scores{$md}[0]\t$pvalue\t$nrand\n"; }

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#       SUBROUTINES
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
sub shuffle (\@) {
	my @ar = @{$_[0]};
	for (my $i = scalar @ar - 1; $i > 0; --$i) {
		my $r = int (rand ($i));
		my $t = $ar[$r];
		$ar[$r] = $ar[$i];
		$ar[$i] = $t; }
	return (@ar); }

