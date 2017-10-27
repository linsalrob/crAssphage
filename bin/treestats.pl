#!/usr/bin/perl -w
use strict;
use Math::Trig;

my $treefile = "";
my $metadata_file = "";
my %metadata_todo = ();
my $min_bs = 0;
my $nrand = 100;
my $weird = 0;
my $extra = 0;
my $usage = "
=== $0 ===

This script calculates clustering statistics of metadata labels in a tree or cladogram.
The statistics are compared to N trees where the leaf labels have been randomized.
The value of the statistic for the unshuffled tree remains the same between multiple runs of the script, the P-value may vary.

(c) Bas E. Dutilh

Use: $0 -t treefile.ph -m metadata.tsv -f field1,field2

    -t   file     File containing tree or cladogram in Newick tree format (option required)
    -m   file     File containing metadata fields per leaf ID in tab-separated format with a header row (option required)
    -c   string   Metadata fields to be tested categorically, multiple separated by commas (at least one of the options -c/-n/-g is required)
    -n   string   Metadata fields to be tested numerically, multiple separated by commas (at least one of the options -c/-n/-g is required)
    -g   string   Metadata fields to be tested as global latitude/longitude coordinates (at least one of the options -c/-n/-g is required)
    -b   float    Minimum branch support value (default: $min_bs)
    -r   integer  Number of randomizations of the leaf labels, determines minimum P-value (default: $nrand)
    -w   binary   Provide metadata file not in tab-separated but in Rob's weird format (default: $weird)
    -x   binary   Extra output for categorical (default: $extra)
    -h            Help: print explanation of the statistics
";

my %explanation = ();
$explanation{"consistency"} = "Average of frequencies of the most frequent metadata annotation in a branch.";
$explanation{"n_perfect"} = "Number of branches where all leaves have the same metadata annotation.";
$explanation{"stand_dev"} = "Average of standard deviations of all metadata annotation values within a branch.";
$explanation{"global_sd"} = "Average of standard deviations of distances between global lat/lon coordinates within a branch.";

ARG: for (my $i = 0; $i < scalar @ARGV; $i += 2) {
        if ($ARGV[$i] eq "-t") {
		if (! (-e $ARGV[$i + 1])) {
			print STDERR "No such file $ARGV[$i + 1] - skipping\n";
			next ARG; }
                $treefile = $ARGV[$i + 1]; }
        elsif ($ARGV[$i] eq "-m") {
                if (! (-e $ARGV[$i + 1])) {
                        print STDERR "No such file $ARGV[$i + 1] - skipping\n";
                        next ARG; }
                $metadata_file = $ARGV[$i + 1]; }
        elsif (($ARGV[$i] eq "-c") || ($ARGV[$i] eq "-n")) {
                foreach my $md (split /,/, $ARGV[$i + 1]) {
                        ++$metadata_todo{$md}{$ARGV[$i]}; } }
	elsif ($ARGV[$i] eq "-g") {
		++$metadata_todo{$ARGV[$i + 1]}{"-g"}; }
	elsif ($ARGV[$i] eq "-b") {
		$min_bs = $ARGV[$i + 1]; }
        elsif ($ARGV[$i] eq "-r") {
		$nrand = $ARGV[$i + 1]; }
	elsif ($ARGV[$i] eq "-w") {
		$weird = $ARGV[$i + 1]; }
	elsif ($ARGV[$i] eq "-x") {
		$extra = $ARGV[$i + 1]; }
	elsif ($ARGV[$i] eq "-h") {
		die "$usage
Categorical statistics:
    consistency   $explanation{\"consistency\"}
    n_perfect     $explanation{\"n_perfect\"}

Numerical statistics:
    stand_dev     $explanation{\"stand_dev\"}

Global coordinates:
    global_sd     $explanation{\"global_sd\"}

"; }
	else {
		print STDERR "Option $ARGV[$i] not recognized, skipping\n"; } }

if ((! (-e $treefile)) || (! (-e $metadata_file)) || (scalar (keys %metadata_todo == 0))) {
	die $usage; }

sub shuffle (\@);
sub average (\@);

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#	READ BRACKETNOTATION FROM TREEFILE
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
my $bracketnotation = "";
open (IN, "<$treefile") or die "Can't open file $treefile: $!\n";
while (my $line = <IN>) {
	chomp ($line);
	$bracketnotation .= $line; }
close (IN);
$bracketnotation =~ s/\s+/ /g;
$bracketnotation =~ s/_R_//g; # THIS LINE IS NEEDED BECAUSE MAFFT ADDS _R_ AT THE START OF THE ID IF SEQUENCE IS REVERSED DURING ALIGNMENT
my $bracketnotation_size = length ($bracketnotation);

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#	DETERMINE ALL PARTITIONS IN BRACKETNOTATION
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
my %partition = ();
for (my $position_rightbracket = 0; $position_rightbracket < $bracketnotation_size; ++$position_rightbracket) {
	if (substr ($bracketnotation, $position_rightbracket, 1) eq ")") {
		FIND_CLOSEST_UNUSED_LEFTBRACKET: for (my $position_leftbracket = $position_rightbracket - 1; $position_leftbracket >= 0; --$position_leftbracket) {
			if (substr ($bracketnotation, $position_leftbracket, 1) eq "(") {
				my $leftbracket_already_used = 0;
				foreach my $used_position_rightbracket (sort keys %partition) {
					if ($partition{$used_position_rightbracket} == $position_leftbracket) {
						++$leftbracket_already_used; } }
				if ($leftbracket_already_used == 0) {
					$partition{$position_rightbracket} = $position_leftbracket;
					last FIND_CLOSEST_UNUSED_LEFTBRACKET; } } } } }
foreach my $position_rightbracket (sort keys %partition) {
        my $bs_and_bl = 1;
        while (substr ($bracketnotation, $position_rightbracket + 1, $bs_and_bl) !~ /[\),;]/) {
                ++$bs_and_bl; }
	$partition{$position_rightbracket} = substr ($bracketnotation, $partition{$position_rightbracket}, $position_rightbracket - $partition{$position_rightbracket} + $bs_and_bl); }

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#	READ METADATA FROM FILE
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
my %metadata = ();
my %all_md_fields = ();
my %cols = ();
open (IN, "<$metadata_file") or die "Can't open file $metadata_file: $!\n$usage";
LINE: while (my $line = <IN>) {
	chomp ($line);
	my @this = split /\t/, $line;

# METADATA IS IN A TAB-DELIMITED TEXT FILE
	if (! $weird) {
		if (scalar keys %cols == 0) {
	    	for (my $v = 0; $v < scalar @this; ++$v) {
	        	$cols{$this[$v]} = $v; }
	        next LINE; }
		else {
			foreach my $v (keys %cols) {
				$metadata{$this[0]}{$v} = $this[$cols{$v}];
				++$all_md_fields{$v}{$this[$cols{$v}]}; } } }

# METADATA IS IN THE OTHER FORMAT
	else {
		my @fields = split / \[/, $this[1];
		for (my $i = 1; $i < scalar @fields; ++$i) {
			$fields[$i] =~ s/\]\s*\Z//;
			my ($k, $v) = split /=/, $fields[$i];
			++$cols{$k};
$this[0] = ">$fields[0]"; # THIS LINE IS NEEDED TO DEAL WITH TREES CONTAINING NAMES FOR LEAVES
			$metadata{$this[0]}{$k} = $v;
			++$all_md_fields{$k}{$v}; } }

	foreach my $md (keys %metadata_todo) {
		if (! (exists ($metadata{$this[0]}{$md}))) {
			die "Error:\tmetadata field '$md' not available for node $this[0] in $metadata_file\n$usage"; } } }

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#	DETERMINE SPECIES AND BRANCH SUPPORT FOR ALL PARTITIONS
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
my %species_in_clades = ();
my %bs = ();
my $max_bs = 0;
foreach my $rb (sort { $a <=> $b; } keys %partition) {
	my $tmp_species = "$partition{$rb}";
	$bs{$rb} = -1;
	if ($tmp_species =~ /.*\)([\-\w\.]+)/) {
		$bs{$rb} = $1; }
	if ($bs{$rb} > $max_bs) {
		$max_bs = $bs{$rb}; }
	$tmp_species =~ s/[\):][\-\w\.]+//g;
	my @species_file = split /[\(\),;]+/, $tmp_species;
	@species_file = sort (@species_file);
	while ($species_file[0] eq "") {
		shift (@species_file); }
	for (my $j = 0; $j < scalar @species_file; ++$j) {
		++$species_in_clades{$rb}{$species_file[$j]}; } }
if (($min_bs != 0) && ($max_bs == 0)) {
	print STDERR "Warning:\tIt appears that $treefile does not contain branch support values -ignoring argument -b ($min_bs)\n"; }
if ($min_bs > $max_bs) {
	die "Error:\tMinimum branch support value ($min_bs) higher than maximum branch support value found in $treefile ($max_bs)\n"; }

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#	LIST ALL SPECIES
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
my @partitions_by_size = sort { length ($partition{$a}) <=> length ($partition{$b}); } keys %partition;
my %species = ();
my @species_list = ();
foreach my $sp (keys %{$species_in_clades{$partitions_by_size[-1]}}) {
	$species{$sp} = $sp;
	if (! (exists ($metadata{$species{$sp}}))) {
		die "Error:\tno metadata for $species{$sp} in $metadata_file\n"; }
	push (@species_list, $sp); }

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#       REMOVE BRANCHES WITH LOW SUPPORT
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
if ($min_bs != 0) {
	my $n_branches = scalar (keys %partition);
	my $n_ignored = 0;
	foreach my $rb (keys %partition) {
		if ($bs{$rb} < $min_bs) {
			delete ($partition{$rb});
			++$n_ignored; } }
	print STDERR "Info:\t$n_ignored/$n_branches branches ingored with support value <$min_bs\n"; }

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#	COUNT METADATA VALUES FOR ALL CLADES
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
my %all_scores = ();
my %best_per_leaf = ();
my %perfect_scores = ();
my %numerical_sd = ();
my %globaldist = ();
# calculate scores nrand times; for the first ($shuf=0) the species are not randomized and for the rest they are
for (my $shuf = 0; $shuf <= $nrand; ++$shuf) {
	# which metadata fields to test?
	foreach my $md (keys %metadata_todo) {
		if (exists ($metadata_todo{$md}{"-c"})) {
			$perfect_scores{$md}[$shuf] = 0;
			# calculate a value for all partitions rb
			foreach my $rb (keys %partition) {
				my %rb_values = ();
				my $n = 0;
				foreach my $sp (keys %{$species_in_clades{$rb}}) {
					++$rb_values{$metadata{$species{$sp}}{$md}};
					++$n; }
				# count frequency of most frequent metadata annotation in this partition
				my @rb_values_sorted = sort { $rb_values{$b} <=> $rb_values{$a}; } keys %rb_values;
				push (@{$all_scores{$md}[$shuf]}, $rb_values{$rb_values_sorted[0]} / $n);
				# extra info: register best score for each species
				foreach my $sp (keys %{$species_in_clades{$rb}}) {
					my $this_score = $rb_values{$metadata{$species{$sp}}{$md}} / $n;
					if ((! exists $best_per_leaf{$md}[$shuf]{$sp}) || ($best_per_leaf{$md}[$shuf]{$sp} < $this_score)) {
						$best_per_leaf{$md}[$shuf]{$sp} = $this_score; } }
				# count if this is a perfect partition
				if ($rb_values{$rb_values_sorted[0]} == $n) {
					++$perfect_scores{$md}[$shuf]; } } }
		elsif (exists ($metadata_todo{$md}{"-n"})) {
			foreach my $rb (keys %partition) {
				# calculate average and standard deviation of values within each partition
				my @values = ();
				foreach my $sp (keys %{$species_in_clades{$rb}}) {
					push (@values, $metadata{$species{$sp}}{$md}); }
				my $this_average = &average (@values);
				my $this_sd = 0;
				foreach my $sp (keys %{$species_in_clades{$rb}}) {
					$this_sd += ($metadata{$species{$sp}}{$md} - $this_average) ** 2; }
				push (@{$numerical_sd{$md}[$shuf]}, sqrt ($this_sd / (scalar (keys %{$species_in_clades{$rb}}) - 1))); } }
		elsif (exists ($metadata_todo{$md}{"-g"})) {
			foreach my $rb (keys %partition) {
				# calculate average and standard deviation of global distances between leaves in a partition
				my @values = ();
				foreach my $sp1 (keys %{$species_in_clades{$rb}}) {
					foreach my $sp2 (keys %{$species_in_clades{$rb}}) {
						my ($lat1, $lon1) = split /,/, $metadata{$species{$sp1}}{$md};
						my ($lat2, $lon2) = split /,/, $metadata{$species{$sp2}}{$md};
						push (@values, &global_distance ($lat1, $lon2, $lat2, $lon2)); } }
				my $this_average = &average (@values);
				my $this_sd = 0;
				foreach my $sp1 (keys %{$species_in_clades{$rb}}) {
					foreach my $sp2 (keys %{$species_in_clades{$rb}}) {
						my ($lat1, $lon1) = split /,/, $metadata{$species{$sp1}}{$md};
						my ($lat2, $lon2) = split /,/, $metadata{$species{$sp2}}{$md};
						$this_sd += (&global_distance ($lat1, $lon2, $lat2, $lon2) - $this_average) ** 2; } }
				push (@{$globaldist{$md}[$shuf]}, sqrt ($this_sd / (scalar (keys %{$species_in_clades{$rb}}) - 1))); } } }

	# after first calculating the values for the unshuffled tree, shuffle the tree for following iterations
	@species_list = shuffle (@species_list);

	my $pos = 0;
	foreach my $sp (keys %species) {
		$species{$sp} = $species_list[$pos];
		++$pos; } }

print "#field\tstatistic\ttype\tvalue\tPval\tpermut\texplanation of the statistic\n";
my $min_pvalue = 1 / $nrand;
foreach my $md (sort keys %metadata_todo) {
	if (exists ($metadata_todo{$md}{"-c"})) {
		my $P_frac = 0;
		my $P_perf = 0;
		my %P_best = ();
		# calculate "mean fraction" statistic of unshuffled tree
		my $frac_unshuf = &average (@{$all_scores{$md}[0]});
		# calculate "average best" statistic of unshuffled tree
		my %md_best_averages = ();
		foreach my $sp (keys %species) {
			$md_best_averages{$metadata{$sp}{$md}}{"t"}[0] += $best_per_leaf{$md}[0]{$sp};
			++$md_best_averages{$metadata{$sp}{$md}}{"n"}[0]; }
		foreach my $field (keys %md_best_averages) {
			$md_best_averages{$field}{"t"}[0] /= $md_best_averages{$field}{"n"}[0]; }
		# calculate P-values
		for (my $shuf = 1; $shuf < $nrand; ++$shuf) {
			# "mean fraction" statistic in shuffled tree
			if (&average (@{$all_scores{$md}[$shuf]}) >= $frac_unshuf) {
				++$P_frac; }
			# "average best" statistic in shuffled tree
			foreach my $sp (keys %species) {
				$md_best_averages{$metadata{$sp}{$md}}{"t"}[$shuf] += $best_per_leaf{$md}[$shuf]{$sp};
				++$md_best_averages{$metadata{$sp}{$md}}{"n"}[$shuf]; }
			foreach my $field (keys %md_best_averages) {
				$md_best_averages{$field}{"t"}[$shuf] /= $md_best_averages{$field}{"n"}[$shuf];
				if ($md_best_averages{$field}{"t"}[$shuf] >= $md_best_averages{$field}{"t"}[0]) {
					++$P_best{$field}; } }
			# "perfect" statistic in shuffled tree
			if ($perfect_scores{$md}[$shuf] >= $perfect_scores{$md}[0]) {
				++$P_perf; } }
		$P_frac /= $nrand;
		if ($P_frac == 0) {
			$P_frac = "<$min_pvalue"; }
		$P_perf /= $nrand;
		if ($P_perf == 0) {
			$P_perf = "<$min_pvalue"; }
		printf "$md\tconsistency\tchar\t%1.5f\t$P_frac\t$nrand\t$explanation{\"consistency\"}\n", $frac_unshuf;
		print "$md\tn_perfect\tchar\t$perfect_scores{$md}[0]\t$P_perf\t$nrand\t$explanation{\"n_perfect\"}\n";

		if ($extra) {
			print "\tmd\tn\tavgbest\tP_best\n";
			foreach my $field (sort keys %P_best) {
				$P_best{$field} /= $nrand;
				if ($P_best{$field} == 0) {
					$P_best{$field} = "<$min_pvalue"; }
				printf "#\t$field\t$all_md_fields{$md}{$field}\t%1.5f\t$P_best{$field}\n", $md_best_averages{$field}{"t"}[0]; } } }

	elsif (exists ($metadata_todo{$md}{"-n"})) {
		my $P_sd = 0;
		# calculate "standard deviation" statistic of unshuffled tree
		my $sd_unshuf = &average (@{$numerical_sd{$md}[0]});
		# calculate P-values
		for (my $shuf = 1; $shuf < $nrand; ++$shuf) {
			# "standard deviation" statistic of shuffled tree
			if (&average (@{$numerical_sd{$md}[$shuf]}) <= $sd_unshuf) {
				++$P_sd; } }
		$P_sd /= $nrand;
		if ($P_sd == 0) {
			$P_sd = "<$min_pvalue"; }
		printf "$md\tstand_dev\tnumer\t%1.2f\t$P_sd\t$nrand\t$explanation{\"stand_dev\"}\n", $sd_unshuf; }

	elsif (exists ($metadata_todo{$md}{"-g"})) {
		my $P_glob = 0;
		# calculate "global_sd" statistic of unshuffled tree
		my $glob_unshuf = &average (@{$globaldist{$md}[0]});
		# calculate P-values
		for (my $shuf = 1; $shuf < $nrand; ++$shuf) {
			# "global_sd" statistic of shuffled tree
			if (&average (@{$globaldist{$md}[$shuf]}) <= $glob_unshuf) {
				++$P_glob; } }
		$P_glob /= $nrand;
		if ($P_glob == 0) {
			$P_glob = "<$min_pvalue"; }
		printf "$md\tglobal_sd\tglobal\t%1.2f\t$P_glob\t$nrand\t$explanation{\"global_sd\"}\n", $glob_unshuf; } }

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

sub average (\@) {
	my $avg = 0;
	foreach my $val (@_) {
		$avg += $val; }
	return ($avg / scalar @_); }

sub global_distance (\@) {
	my $phi1 = deg2rad (90 - $_[0]); my $theta1 = deg2rad ($_[1]);
	my $phi2 = deg2rad (90 - $_[2]); my $theta2 = deg2rad ($_[3]);
	if (($phi1 == $phi2) && ($theta1 == $theta2)) {
		return (0); }
	return (6373 * acos (sin ($phi1) * sin ($phi2) * cos ($theta1 - $theta2) + cos ($phi1) * cos ($phi2))); }

