#!/bin/bash

###############################################################################################################
#                                                                                                             #
# Run gretel on a bunch of files. This assumes we've already run bowtie2 and have indexed bam files.          #
#                                                                                                             #
# This is just a wrapper around GRETEL code from Sam Nicholls, and is specific for crAssphage.                #
#                                                                                                             #
# Written by Rob Edwards. (c) Rob Edwards 2015-2018. Use at your own risk!                                    #
#                                                                                                             #
#                                                                                                             #
###############################################################################################################

# list of SRR IDs that we need to analyze (just the ID)
LISTOFFILES=bowtie.txt
# where to find the bam file
BAMDIR=bowtie
# where to put the output
GRETELDIR=GRETEL
# where is crAssphage
CRASSPHAGE=$HOME/phage/crAssphage/JQ995537.fna

export PATH=$PATH:$HOME/opt/bin/
SRR=$(head -n $SGE_TASK_ID $LISTOFFILES | tail -n 1);


if [ -e $GRETELDIR/$SRR ]; then
	echo "$GRETELDIR/$SRR already exists. Exiting"
	exit $E_BADARGS
fi

mkdir -p $GRETELDIR/$SRR
cd $GRETELDIR/$SRR
ln -s ../../$BAMDIR/$SRR.bam .

BAM=$SRR.bam
if [ ! -h $BAM ]; then
	echo "NO BAM FILE FOR '$BAM' in $PWD";
	exit $E_BADARGS
fi

if [ ! -e ../../$BAMDIR/$SRR.bam.bai ]; then samtools index ../../$BAMDIR/$SRR.bam; fi
ln -s ../../$BAMDIR/$SRR.bam.bai .


python2.7 $HOME/GitHubs/gretel-test/snpper.py -b $BAM -r JQ995537 -s 25633 -e 26964 > ${SRR}.A.vcf
python2.7 $HOME/GitHubs/gretel-test/snpper.py -b $BAM -r JQ995537 -s 33708 -e 35062 > ${SRR}.B.vcf
python2.7 $HOME/GitHubs/gretel-test/snpper.py -b $BAM -r JQ995537 -s 43819 -e 45057 > ${SRR}.C.vcf

for P in A B C; do 
	bgzip ${SRR}.$P.vcf
	tabix ${SRR}.$P.vcf.gz
done



mkdir gretel

gretel --master $CRASSPHAGE -s 25633 -e 26964 -p 10000 -o gretel $BAM ${SRR}.A.vcf.gz  JQ995537  > gretel/gretel_A_output 2> gretel/gretel_A_err
if [ -e gretel/out.fasta ]; then mv gretel/out.fasta gretel/pcrA.fasta; fi
mv gretel/gretel.crumbs gretel/gretel.A.crumbs
gretel --master $CRASSPHAGE -s 33708 -e 35062 -p 10000 -o gretel $BAM ${SRR}.B.vcf.gz  JQ995537  > gretel/gretel_B_output 2> gretel/gretel_B_err
if [ -e gretel/out.fasta ]; then mv gretel/out.fasta gretel/pcrB.fasta; fi
mv gretel/gretel.crumbs gretel/gretel.B.crumbs
gretel --master $CRASSPHAGE -s 43819 -e 45057 -p 10000 -o gretel $BAM ${SRR}.C.vcf.gz  JQ995537  > gretel/gretel_C_output 2> gretel/gretel_C_err
if [ -e gretel/out.fasta ]; then mv gretel/out.fasta gretel/pcrC.fasta; fi
mv gretel/gretel.crumbs gretel/gretel.C.crumbs

