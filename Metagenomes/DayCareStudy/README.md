# Day Care Study

We are grateful to Aaron J. Prussin II, Pedro J. Torres, John Shimashita, Steven R. Head, Scott T. Kelley, and Linsey C. Marr for allowing us to include this data prior to publication of their paper.

This data comes from their study looking at the microbiome of the built environment, specifically a day care nursery that has a lot of crAssphage flying around!

# Data processing

Raw fastq files were cleaned using [prinseq++](https://github.com/Adrian-Cantu/PRINSEQ-plus-plus) with the following command:

```
prinseq++ -threads 20 -fastq raw_reads/AP-DNA-16_S8_R1_001.fastq.gz -fastq2 raw_reads/AP-DNA-16_S8_R2_001.fastq.gz -out_name qc/AP-DNA-16_S8  -min_qual_mean 20 -ns_max_n 0 -derep -trim_qual_right=20 -lc_entropy -min_len 30
```

They were aligned to crAssphage with bowtie2 using the following command:

```
bowtie2 -p 8 -q --no-unal --no-unal -x /home3/redwards/Phage/crAssphage/JQ995537 -1 qc/AP-DNA-16_S8_good_out_R1.fastq -2 qc/AP-DNA-16_S8_good_out_R2.fastq -U qc/AP-DNA-16_S8_single_out_R1.fastq -U qc/AP-DNA-16_S8_single_out_R2.fastq  | samtools view -bS - | samtools sort -o rebowtie/$F.bam
```

Next, haplotypes were called using gretel:

```
# first we make fake SNP variants and compress and index them
python2.7 $HOME/GitHubs/gretel-test/snpper.py -b $BAM -r JQ995537 -s 25633 -e 26964 > ${SRR}.A.vcf
python2.7 $HOME/GitHubs/gretel-test/snpper.py -b $BAM -r JQ995537 -s 33708 -e 35062 > ${SRR}.B.vcf
python2.7 $HOME/GitHubs/gretel-test/snpper.py -b $BAM -r JQ995537 -s 43819 -e 45057 > ${SRR}.C.vcf

for P in A B C; do 
	bgzip ${SRR}.$P.vcf
	tabix ${SRR}.$P.vcf.gz
done
```

Next, we use gretel to call the variants.


```
gretel --master $CRASSPHAGE -s 25633 -e 26964 -p 10000 -o gretel $BAM ${SRR}.A.vcf.gz  JQ995537  > gretel/gretel_A_output 2> gretel/gretel_A_err
if [ -e gretel/out.fasta ]; then mv gretel/out.fasta gretel/pcrA.fasta; fi
mv gretel/gretel.crumbs gretel/gretel.A.crumbs
gretel --master $CRASSPHAGE -s 33708 -e 35062 -p 10000 -o gretel $BAM ${SRR}.B.vcf.gz  JQ995537  > gretel/gretel_B_output 2> gretel/gretel_B_err
if [ -e gretel/out.fasta ]; then mv gretel/out.fasta gretel/pcrB.fasta; fi
mv gretel/gretel.crumbs gretel/gretel.B.crumbs
gretel --master $CRASSPHAGE -s 43819 -e 45057 -p 10000 -o gretel $BAM ${SRR}.C.vcf.gz  JQ995537  > gretel/gretel_C_output 2> gretel/gretel_C_err
if [ -e gretel/out.fasta ]; then mv gretel/out.fasta gretel/pcrC.fasta; fi
mv gretel/gretel.crumbs gretel/gretel.C.crumbs
```

This results in the three files:

# gretel_strains_pcrA.fasta
# gretel_strains_pcrB.fasta
# gretel_strains_pcrC.fasta



