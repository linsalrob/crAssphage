#!/bin/bash
source /etc/profile.d/modules.sh
module load samtools
module load bcftools
for f in *.sam
do
	samtools view -bS $f > ${f%%.*}.bam
	echo "Changed $f to BAM"	
	samtools sort ${f%%.*}.bam -o ${f%%.*}.sorted.bam
	echo "Sorted ${f%%.*} BAM"
	samtools index ${f%%.*}.sorted.bam #this produces S1CA_phi6.sorted.bai
	echo "Indexed ${f%%.*} BAM"
done
