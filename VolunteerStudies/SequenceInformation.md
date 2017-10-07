# Sequence Information

We have generated a lot of sequences from individuals, and we have sequenced both ends of the reads. This page just 
keeps track of the sequences and how we have processed them.

The table below shows the current status of the sequences. For each sample, we try and sequence from 
the left and right ends. We then manually trim those using FinchTV, and align them using mafft with this command:

```angular2html
export SID=DE12
mafft --adjustdirectionaccurately  --maxiterate 16   --thread 8 $SID.fna > $SID.aln; clustalx $SID.aln
clustalx $SID.aln
perl ../bin/fasta_consensus.pl $SID.aln  > ../final_sequences/$SID.fna
```

The `clustalx` step in there is a manual inspection of the alignment to ensure that we have included enough
 sequence from each read to generate a meaningful overlap.
 
The data in the table shows the sample ID, the identifiers of the left and right reads, the volunteer ID and sample 
number, the type of sample (whether it was a pre-enrollment, daily, or weekly sample), the date the sample was collected
or processed, and whether the sequencing
was successful or not (this is basically from our manual inspection of the ab1 files), and we have left and
right reads respectively. (Note: We are not sure why the right side appears to fail a lot more than the left side. 
It suggests there is an issue with that primer.) If we have a failed read, we can't generate a consensus sequence
and only have a single read. These are placed in the [Singles](Sequences/raw_data/final_sequences/Singles) directory.
If we have left and right reads, we grab the consensus sequence, (using [fasta_consensus.pl](../bin/fasta_consensus.pl)),
and then we have the final sequence length.

Note that the ab1 files are in the two tarballs [ab1.tar.bz2](Sequences/raw_data/ab1.tar.bz2) and 
[Failed.tar.bz2](Sequences/raw_data/Failed.tar.bz2). The trimmed fasta files are in the tarball 
[fasta_trimmed]

Label | Left read | Right read | Volunteer ID | Sample | Sample Type | Date collected |  Status  |  Final Status  |  Length (bp) 
--- | --- | --- | --- | --- | --- | --- | --- | --- | ---
D1V11 | UB_019_D1V11F_V11F_A04.ab1 | UB_020_D1V11R_V11R_B04.ab1 | D | 1 | Pre | 20170421 |  Failed / Failed  |  -  |  -
DD10V12 |  |  | D | 8 | Daily  | 20170707 |  |  | 
DD11V12 |  |  | D | 9 | Daily  | 20170710 |  |  | 
DD3V11 | UB_015_DD3V11F_V11F_E03.ab1 | UB_016_DD3V11R_V11R_F03.ab1 | D | 3 | Daily | 20170622 |  Passed / Failed  |  Single read  | 710
DD5V11 | UB_001_DD5V11F_V11F_G01.ab1 | UB_002_DD5V11R_V11R_H01.ab1 | D | 4 | Daily | 20170623 |  Passed / Failed  |  Single read  | 421
DD7V12 |  |  | D | 5 | Daily  | 20170704 |  |  | 
DD8V12 |  |  | D | 6 | Daily  | 20170705 |  |  | 
DD9V12 |  |  | D | 7 | Daily  | 20170706 |  |  | 
DE10V12 |  |  | E | 8 | Daily  | 20170714 |  |  | 
DE11V12 | DA_039_DE11V12F_V12F_G05.ab1 | DA_040_DE11V12R_V12R_H05.ab1 | E | 9 | Daily  | 20170717 |  Passed / Failed  |  Single read  | 890
DE12V12 | DA_033_DE12V12F_V12F_A05.ab1 | DA_034_DE12V12R_V12R_B05.ab1 | E | 10 | Daily  | 20170718 |  Passed / Passed  |  Consensus  | 1123
DE13V12 | DA_021_DE13V12F_V12F_E03.ab1 | DA_022_DE13V12R_V12R_F03.ab1 | E | 11 | Daily  | 20170719 |  Passed / Passed  |  Consensus  | 1135
DE14V12 | DA_035_DE14V12F_V12F_C05.ab1 | DA_036_DE14V12R_V12R_D05.ab1 | E | 12 | Daily  | 20170720 |  Passed / Passed  |  Consensus  | 1135
DE15V12 | DA_029_DE15V12F_V12F_E04.ab1 | DA_030_DE15V12R_V12R_F04.ab1 | E | 13 | Daily  | 20170721 |  Passed / Passed  |  Consensus  | 1133
DE6V12 |  |  | E | 4 | Daily  | 20170710 |  |  | 
DE7V12 |  |  | E | 5 | Daily  | 20170711 |  |  | 
DE8V12 |  |  | E | 6 | Daily  | 20170712 |  |  | 
DE9V12 |  |  | E | 7 | Daily  | 20170713 |  |  | 
E1V11 | UB_017_E1V11F_V11F_G03.ab1 | UB_018_E1V11R_V11R_H03.ab1 | E | 1 | Pre | 20170421 |  Passed / Passed  |  Consensus  | 1241
Ke1V11 | UB_005_Ke1V11F_V11F_C02.ab1 | UB_006_Ke1V11R_V11R_D02.ab1 | Ke | 1 | Pre | 20170525 |  Passed / Failed  |  Single read  | 487
MiV11 | UB_011_MiV11F_V11F_A03.ab1 | UB_012_MiV11R_V11R_B03.ab1 | Mi | 1 | Pre | 20170524 |  Passed / Failed  |  Single read  | 470
WD1V11 | UB_021_WD1V11F_V11F_C04.ab1 | UB_022_WD1V11R_V11R_D04.ab1 | D | 1 | Weekly | 20170601 |  Passed /Passed  |   Consensus  | 1242
WD2V11 | UB_009_WD2V11F_V11F_G02.ab1 | UB_010_WD2V11R_V11R_H02.ab1 | D | 2 | Weekly | 20170608 |  Passed / Passed  |  Consensus  | 1285
WD6V12 | DA_027_WD6V12F_V12F_C04.ab1 | DA_028_WD6V12R_V12R_D04.ab1 | D | 10 | Weekly | 20170722 |  Passed / Failed  |  Single Read  | 780
WD7V12 | DA_031_WD7V12F_V12F_G04.ab1 | DA_032_WD7V12R_V12R_H04.ab1 | D | 11 | Weekly | 20170729 |  Passed / Passed  |  Consensus  | 1134
WE1V11 | UB_023_WE1V11F_V11F_E04.ab1 | UB_024_WE1V11R_V11R_F04.ab1 | E | 1 | Weekly | 20170522 |  Passed / Passed  |  Consensus  | 1269
WE2V11 | UB_007_WE2V11F_V11F_E02.ab1 | UB_008_WE2V11R_V11R_F02.ab1 | E | 2 | Weekly | 20170601 |  Passed / Passed  |  Consensus  | 1267
WE5V12 | DA_023_WE5V12F_V12F_G03.ab1 | DA_024_WE5V12R_V12R_H03.ab1 | E | 3 | Weekly | 20170629 |  Passed / Passed  |  Consensus  | 1139
WE6V12 | DA_025_WE6V12F_V12F_A04.ab1 | DA_026_WE6V12R_V12R_B04.ab1 | E | 15 | Weekly | 20170731 |  Passed / Passed  |  Consensus  | 1134
WE7V12 | DA_037_WE7V12F_V12F_E05.ab1 | DA_038_WE7V12R_V12R_F05.ab1 | E | 14 | Weekly | 20170725 |  Passed / Passed  |  Consensus  | 1145
WMi4V11 | UB_003_WMi4V11F_V11F_A02.ab1 | UB_004_WMi4V11R_V11R_B02.ab1 | Mi | 4 | Weekly | 20170517 |  Failed / Failed  |  -  |  -
WMiV11 | UB_013_WMiV11F_V11F_C03.ab1 | UB_014_WMiV11R_V11R_D03.ab1 | Mi | 3 | Weekly | 20170622 |  Failed / Failed  |  -  |  -
Y1V11 | UB_025_Y1V11F_V11F_G04.ab1 | UB_026_Y1V11R_V11R_H04.ab1 | Y | 1 | Pre | 20170518 |  Passed / Passed  |  Consensus  | 1257