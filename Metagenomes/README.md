# crAssphage sequences extracted from metagenomes

We extracted all the crAssphage sequences from all the metagenomes in the [sequence read archive](https://www.ncbi.nlm.nih.gov/sra/), and have built trees based on those sequences.

We have a complex pipeline, and in other work Rob is collaborating with researchers at Indiana University to make this a widely applicable pipeline that you can use for any sequence. However, at the moment we don't have a working pipeline that anyone can use. 

Here are the steps that we use, and the commands that were performed to.

## Step 1: Map all SRA Runs to crAssphage

We have developed a protocol for separating amplicon, metagenome, and other data sets from the SRA, called [PARTIE](https://github.com/linsalrob/partie). We start with this list of WGS metagenomes, and compare 100,000 reads from each metagenome to crAssphage using bowtie2. Our bowtie command is:

```
bowtie2 -p 6 -q --no-unal -x $1 -U $SRR | samtools view -bS - | samtools sort - output_dir/$SRR
samtools index output_dir/$SRR
```

(Note: This is probably incorrect, and we should adjust these as described [by Sam Nicholls](https://samnicholls.net/2016/12/24/bowtie2-metagenomes/))

Once we have a set of metagenomes where _any_ read in the 100,000 reads matches crAssphage, we go back and remap the whole run to crAssphage. This is a heuristic that allows us to search the whole SRA in a reasonable time frame.

At the end of this, we have a series of .bam files, one per metagenome.

## Step 2. Pull out those sequences that map to the three PCR regions. 

We use either the code we [wrote in python](https://github.com/linsalrob/EdwardsLab/blob/master/crAssphage/extract_pcr_regions.py) or this samtools command to pull out the reads from the bam file that overlap our PCR regions (bamfiles is a directory with all the bamfiles).

```
for BAM in $(ls bamfiles | grep bam$):
do 
	FQ=$(echo $BAM|sed -e 's/bam/fastq/');
	samtools view ../../bamfiles/$BAM JQ995537:25633-26964 -b | samtools fastq - > fastq_A/$FQ
	samtools view ../../bamfiles/$BAM JQ995537:33708-35062 -b | samtools fastq - > fastq_B/$FQ
	samtools view ../../bamfiles/$BAM JQ995537:43819-45057 -b | samtools fastq - > fastq_C/$FQ
done
```

This creates the fastq files in fastq_A, fastq_B, and fastq_C that are provided here as bzip2 tar files.

## Step 3. Assemble those sequences

We compared using cap3, spades using -k77,99 or -k21,33, and newbler with the commands below to assemble the sequences. Note that cap3 and newbler used fasta, spades used fastq.

```
cap3 fasta_A/DRR002659.fasta -p 90 -o 20 -k 0 -h 70 -s 251 -i 30 -j 31
spades.py -k 77,99 -s fastq_A/DRR002659.fastq -o DRR002659/
spades.py -k 21,33 -s fastq_A/DRR002659.fastq -o DRR002659/
runAssembly -o DRR002659/ -p fasta_A/DRR002659.fasta
```

We found that cap3 gave lots of short contigs, spades kept crashing with core dumps, and newbler gave a reasonable number of long contigs, so we used the newbler data.

## Step 4. Extract the country information from the SRA metadata.

We downloaded the [biosample xml file from NCBI](ftp://ftp.ncbi.nih.gov/biosample) and then [parsed the XML using perl](../bin/parse_sra_metadata.pl). (Note that we don't read the whole file in because it is about 15G, so we just parse line by line.)

We also get the SRR data from the SQL using this code:

```
for SRR in $(cat /home3/redwards/Phage/crAssphage/SRA_alignments/crAssphage_srr_ids.txt);
do 
	SRS=$(sqlite3 SRAmetadb.sqlite "select sample_accession, xref_link from sample where sample_accession in (select sample_accession from experiment where experiment_accession in (select experiment_accession from run where run_accession='$SRR'))");
	echo -e "$SRR\t$SRS";
done > /home3/redwards/Phage/crAssphage/SRA_alignments/sample_ids.txt
```
(Note, you can learn more about the SRA metadata from [our blogposts on the SRA](https://edwards.sdsu.edu/SRA))

Then we just join those two tables to get the SRA information we need.

## Step 5. Align the assemblies

We renamed the assemblies, and added the country information for each of the reads from step 4, and then put this through our usual pipeline with a couple of minor modifications:
1. We switched the phylogenetics program from PHYLIP to FastTreeDbl to handle the larger datasets
2. We use the heuristic tree parsing code (rename_tree_fast.py) that uses regular expressions to rename the tree rather than parsing the whole tree.

That approach gives the trees in the repostiory.




