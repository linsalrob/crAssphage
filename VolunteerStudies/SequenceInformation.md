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
 
The data in the table shows the sample ID, the identifiers of the left and right reads, whether the sequencing
was successful or not (this is basically from our manual inspection of the ab1 files), and we have left and
right reads respectively. (Note: We are not sure why the right side appears to fail a lot more than the left side. 
It suggests there is an issue with that primer.) If we have a failed read, we can't generate a consensus sequence
and only have a single read. These are placed in the [Singles](Sequences/raw_data/final_sequences/Singles) directory.
If we have left and right reads, we grab the consensus sequence, (using [fasta_consensus.pl](../bin/fasta_consensus.pl)),
and then we have the final sequence length.

Note that the ab1 files are in the two tarballs [ab1.tar.bz2](Sequences/raw_data/ab1.tar.bz2) and 
[Failed.tar.bz2](Sequences/raw_data/Failed.tar.bz2). The trimmed fasta files are in the tarball 
[fasta_trimmed]

Sample | Left ID | Right ID | Status | Final Status | Length (bp) |
--- | --- | --- | --- | --- | --- |
D1 | UB_019_D1V11F_V11F_A04 | UB_020_D1V11R_V11R_B04 | Failed / Failed | - | -
DD3 | UB_015_DD3V11F_V11F_E03 | UB_016_DD3V11R_V11R_F03 | Passed / Failed | Single read | 710
DD5 | UB_001_DD5V11F_V11F_G01 | UB_002_DD5V11R_V11R_H01 | Passed / Failed | Single read | 421
DE11 | DA_039_DE11V12F_V12F_G05 | DA_040_DE11V12R_V12R_H05 | Passed / Failed | Single read | 890
DE12 | DA_033_DE12V12F_V12F_A05 | DA_034_DE12V12R_V12R_B05 | Passed / Passed | Consensus | 1123
DE13 | DA_021_DE13V12F_V12F_E03 | DA_022_DE13V12R_V12R_F03 | Passed / Passed | Consensus | 1135 
DE14 | DA_035_DE14V12F_V12F_C05 | DA_036_DE14V12R_V12R_D05 | Passed / Passed | Consensus | 1135
DE15 | DA_029_DE15V12F_V12F_E04 | DA_030_DE15V12R_V12R_F04 | Passed / Passed | Consensus | 1133 
E1 | UB_017_E1V11F_V11F_G03 | UB_018_E1V11R_V11R_H03 | Passed / Passed | Consensus | 1241
Ke1 | UB_005_Ke1V11F_V11F_C02 | UB_006_Ke1V11R_V11R_D02 | Passed / Failed | Single read | 487
Mi | UB_011_MiV11F_V11F_A03 | UB_012_MiV11R_V11R_B03 | Passed / Failed | Single read | 470
WD1 | UB_021_WD1V11F_V11F_C04 | UB_022_WD1V11R_V11R_D04 | Passed /Passed |  Consensus | 1242
WD2 | UB_009_WD2V11F_V11F_G02 | UB_010_WD2V11R_V11R_H02 | Passed / Passed | Consensus | 1285
WD6 | DA_027_WD6V12F_V12F_C04 | DA_028_WD6V12R_V12R_D04 | Passed / Failed | Single Read | 780
WD7 | DA_031_WD7V12F_V12F_G04 | DA_032_WD7V12R_V12R_H04 | Passed / Passed | Consensus | 1134
WE1 | UB_023_WE1V11F_V11F_E04 | UB_024_WE1V11R_V11R_F04 | Passed / Passed | Consensus | 1269
WE2 | UB_007_WE2V11F_V11F_E02 | UB_008_WE2V11R_V11R_F02 | Passed / Passed | Consensus | 1267
WE5 | DA_023_WE5V12F_V12F_G03 | DA_024_WE5V12R_V12R_H03 | Passed / Passed | Consensus | 1139
WE6 | DA_025_WE6V12F_V12F_A04 | DA_026_WE6V12R_V12R_B04 | Passed / Passed | Consensus | 1134
WE7 | DA_037_WE7V12F_V12F_E05 | DA_038_WE7V12R_V12R_F05 | Passed / Passed | Consensus | 1145
WMi4 | UB_003_WMi4V11F_V11F_A02 | UB_004_WMi4V11R_V11R_B02 | Failed / Failed | - | -
WMi | UB_013_WMiV11F_V11F_C03 | UB_014_WMiV11R_V11R_D03 | Failed / Failed | - | -
Y1 | UB_025_Y1V11F_V11F_G04 | UB_026_Y1V11R_V11R_H04 | Passed / Passed | Consensus | 1257

