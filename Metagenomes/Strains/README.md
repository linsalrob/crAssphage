# Strain-resolved bioinformatics

For each of the metagenomes where we have identified crAssphage, we attempted to generate the strains ('haplotypes') within each sequence using [gretel](https://gretel.readthedocs.io) by Sam Nicholls. ([Here is the borxiv paper](http://www.biorxiv.org/content/early/2016/08/02/067215) describing gretel.)

We started by aligning all the reads from the SRA library against crAssphage, and then used snpper.py from [gretel-test](https://github.com/SamStudio8/gretel-test) to create VCF (variant files) for each of the three regions that we are interested in:

In this case, SRR is our SRA run identifier, and of course we do this in a loop (or using our cluster).

```
python2.7 snpper.py -b ${SRR}.bam -r JQ995537 -s 25633 -e 26964 > ${SRR}.A.vcf
python2.7 snpper.py -b ${SRR}.bam -r JQ995537 -s 33708 -e 35062 > ${SRR}.B.vcf
python2.7 snpper.py -b ${SRR}.bam -r JQ995537 -s 43819 -e 45057 > ${SRR}.C.vcf
```

Then, for each of these vcf files, we need to bgzip them and index them:

```
for P in A B C; do
	bgzip ${SRR}.$P.vcf
	tabix ${SRR}.$P.vcf.gz
done
```

And then we can run gretel on each of our variant calling files like this:

```
gretel --master JQ995537.fna -s 25633 -e 26964 -o gretel $BAM ${SRR}.A.vcf.gz  JQ995537  > gretel/gretel_A_output 2> gretel/gretel_A_err
if [ -e gretel/out.fasta ]; then mv gretel/out.fasta gretel/pcrA.fasta; fi
mv gretel/gretel.crumbs gretel/gretel.A.crumbs
gretel --master JQ995537.fna -s 33708 -e 35062 -o gretel $BAM ${SRR}.B.vcf.gz  JQ995537  > gretel/gretel_B_output 2> gretel/gretel_B_err
if [ -e gretel/out.fasta ]; then mv gretel/out.fasta gretel/pcrB.fasta; fi
mv gretel/gretel.crumbs gretel/gretel.B.crumbs
gretel --master JQ995537.fna -s 43819 -e 45057 -o gretel $BAM ${SRR}.C.vcf.gz  JQ995537  > gretel/gretel_C_output 2> gretel/gretel_C_err
if [ -e gretel/out.fasta ]; then mv gretel/out.fasta gretel/pcrC.fasta; fi
mv gretel/gretel.crumbs gretel/gretel.C.crumbs
```

Finally, we combine all those outputs into a single fasta file for each of the three primers, that we've called gretel_strains_pcrA.fasta, gretel_strains_pcrB.fasta, and gretel_strains_pcrC.fasta and you can access in [PrimerA](PrimerA/), [PrimerB](PrimerB/), and [PrimerC](PrimerC/).

In those directories, as usual, we have Makefiles for aligning the sequences and building trees. 





