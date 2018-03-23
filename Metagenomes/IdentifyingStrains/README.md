# Identifying strains from metagenomes

One of our key approaches to analyze crAssphage from around the globe is to identify strains from metagenomes. There are several tools that purport to be able to do this, including [DESMAN](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1309-9) and [CONSTRAINS](https://www.nature.com/articles/nbt.3319), but at the time we started this analysis the _only_ working tool that gave us reproducible results was GRETEL whose papers are on the biorxiv: [Probabilistic Recovery Of Cryptic Haplotypes From Metagenomic Data](https://www.biorxiv.org/content/early/2017/03/17/117838) and [Computational haplotype recovery and long-read validation identifies novel isoforms of industrially relevant enzymes from natural microbial communities](https://www.biorxiv.org/content/early/2018/01/13/223404). 

Gretel allowed us to specify the three regions that we want to recover strains from, allowed us to identify strains from bam files, and allowed us to compute this on our cluster. Note thatthis does not mean anything about the other tools, just that when we started this analysis (in 2015/2016) they were still in development and we could not get them to work. I know that in the intervening years they (and others) have been developed further. Therefore, our pipeline uses GRETEL for all the analysis. 

Our code is also built around having a cluster with [Sun Grid Engine](https://en.wikipedia.org/wiki/Oracle_Grid_Engine) for job submission. The key piece that you will see in our code is the environment variable `$SGE_TASK_ID` which is a unique number that each job gets. This allows us to process hundreds or thousands of sequences with a single command. Essentially, `$SGE_TASK_ID` becomes a single integer from 1 to however many data sets we are analyzing.

## Quality Control

We start by using [prinseq++](https://github.com/Adrian-Cantu/PRINSEQ-plus-plus) to perform quality contol on the sequences. We trim the sequences using these parameters:

```
prinseq++ -threads 20 -fastq R1.fastq -fastq2 R2.fastq -out_name trimmed/data -DNA-16_S8  -min_qual_mean 20 -ns_max_n 0 -derep -trim_qual_right=20 -lc_entropy -min_len 30
```

## Mapping the reads

Next, we map the reads using bowtie2. We typically use the code in [bowtie.sh](bowtie.sh) to map reads. The complexity in that code is figuring out if we have paired end reads, singletons (usually from the cleaning step), or unpaired reads. The essential bowtie2 command we use is:

```
bowtie2 -p 8 -q --no-unal --no-unal -x /home3/redwards/Phage/crAssphage/JQ995537 -1 trimmed/R1.fastq -2 trimmed/R2.fastq | samtools view -bS - | samtools sort -o bam/FILE.bam
```

### Running Gretel

Next, we run [gretel.sh](gretel.sh) to run gretel for each of our three primer regions that we are interested in.

We start by creating [VCF](http://www.internationalgenome.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-40/) files for each of the three regions based on the bam files. Note that the gretel-test code was provided by Sam Nicholls.

In this code `$BAM` is our BAM file, and `$SRR` is the ID of the file that we want to identify variants in.

```
python2.7 $HOME/GitHubs/gretel-test/snpper.py -b $BAM -r JQ995537 -s 25633 -e 26964 > ${SRR}.A.vcf
python2.7 $HOME/GitHubs/gretel-test/snpper.py -b $BAM -r JQ995537 -s 33708 -e 35062 > ${SRR}.B.vcf
python2.7 $HOME/GitHubs/gretel-test/snpper.py -b $BAM -r JQ995537 -s 43819 -e 45057 > ${SRR}.C.vcf
```

We then compress and index those files

```
for P in A B C; do 
	bgzip ${SRR}.$P.vcf
	tabix ${SRR}.$P.vcf.gz
done
```

And we run gretel on each of the three variants, being sure to save the crumbs and results files:

```
gretel --master $CRASSPHAGE -s 25633 -e 26964 -p 10000 -o gretel $BAM ${SRR}.A.vcf.gz  JQ995537  > gretel/gretel_A_output 2> gretel/gretel_A_err
if [ -e gretel/out.fasta ]; then mv gretel/out.fasta gretel/pcrA.fasta; fi
mv gretel/gretel.crumbs gretel/gretel.A.crumbs
```


Now we have the variants in regions we call _pcrA.fasta_ (and of course, we generate files for pcrB and pcrC).


