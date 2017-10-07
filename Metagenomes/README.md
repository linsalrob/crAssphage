# crAssphage sequences extracted from metagenomes

We extracted all the crAssphage sequences from all the metagenomes in the [sequence read archive](https://www.ncbi.nlm.nih.gov/sra/), and have built trees based on those sequences.

We have a complex pipeline, and in other work Rob is collaborating with researchers at Indiana University to make this a widely applicable pipeline that you can use for any sequence. However, at the moment we don't have a working pipeline that anyone can use. 

Here are the steps that we use, and the commands that were performed to.

## Map all SRA Runs to crAssphage

We have developed a protocol for separating amplicon, metagenome, and other data sets from the SRA, called [PARTIE](https://github.com/linsalrob/partie). We start with this list of WGS metagenomes, and compare 100,000 reads from each metagenome to crAssphage using bowtie2. Our bowtie command is:

```
bowtie2 -p 6 -q --no-unal -x $1 -U $SRR | samtools view -bS - | samtools sort - output_dir/$SRR
samtools index output_dir/$SRR
```

(Note: This is probably incorrect, and we should adjust these as described [by Sam Nicholls](https://samnicholls.net/2016/12/24/bowtie2-metagenomes/))

Once we have a set of metagenomes where _any_ read in the 100,000 reads matches crAssphage, we go back and remap the whole run to crAssphage. This is a heuristic that allows us to search the whole SRA in a reasonable time frame.

At the end of this, we have a series of .bam files, one per metagenome. Please note, we do not include these bam files in the repository because they are too large.

These files are used to extract the [haplotypes](Haplotypes/) using Gretel.


## Attach the metadata for the metagenomes

We downloaded the [biosample xml file from NCBI](ftp://ftp.ncbi.nih.gov/biosample) and then [parsed the XML using perl](../bin/parse_sra_metadata.pl). (Note that we don't read the whole file in because it is about 15G, so we just parse line by line.)

We also get the SRR data from the SQL using this code (crAssphage_srr_ids.txt is just a file that has all the SRR ids for the runs that match to crAssphage, one per line). SRAmetadb.sqlite is the [SRA metadata](https://edwards.sdsu.edu/research/sra-metadata/) database.

```
for SRR in $(cat crAssphage_srr_ids.txt);
do 
	SRS=$(sqlite3 SRAmetadb.sqlite "select sample_accession, xref_link from sample where sample_accession in (select sample_accession from experiment where experiment_accession in (select experiment_accession from run where run_accession='$SRR'))");
	echo -e "$SRR\t$SRS";
done > /home3/redwards/Phage/crAssphage/SRA_alignments/sample_ids.txt
```
(Note, you can learn more about the SRA metadata from [our blogposts on the SRA](https://edwards.sdsu.edu/SRA))

Then we just join those two tables to get the SRA information we need.


## Other Metagenome Sequences

We have also extracted sequences from other metagenomes. 


