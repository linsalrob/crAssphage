
########################################################################################
#                                                                                      #
# Run bowtie 2 on a bunch of sequence files where we may have left, right, or single   #
# reads. We figure out which files we have, and then call the bowtie2 command.         #
#                                                                                      #
# This is designed to by run on a cluster using SGE.                                   #
#                                                                                      #
# (c) Rob Edwards, 2015-2018. Written by Rob Edwards                                   #
#                                                                                      #
########################################################################################


# list of SRA IDs to compute on
SRRFILE=seqids.txt
# directory with those files
CT=prinseqpp_output

# output directory
OD=bowtie
mkdir -p $OD

F=$(head -n $SGE_TASK_ID $SRRFILE| tail -n 1)
echo "File is $F";
echo "Job is $SGE_TASK_ID"


if [ -e $OD/$F.bam ]; then
	echo "$OD/$F.bam was already found. Skipped"
	exit $E_BADARGS
fi


export PATH=$PATH:$HOME/opt/bin/:$HOME/opt/bowtie2/current/

if [ -e $CT/${F}_good_out_R1.fastq ] && [ -e $CT/${F}_good_out_R2.fastq ]; then
	O=" -1 $CT/${F}_good_out_R1.fastq -2 $CT/${F}_good_out_R2.fastq";
elif [ -e $CT/${F}_good_out_R1.fastq ]; then
	O=" -U $CT/${F}_good_out_R1.fastq";
elif [ -e $CT/${F}_good_out_R2.fastq ]; then
	O=" -U $CT/${F}_good_out_R2.fastq";
fi

if [ -e $CT/${F}_good_out.fastq ]; then
	O=$O" -U $CT/${F}_good_out.fastq";
fi


if [ -e $CT/${F}_single_out_R1.fastq ]; then
	O=$O" -U $CT/${F}_single_out_R1.fastq";
fi

if [ -e $CT/${F}_single_out_R2.fastq ]; then
	O=$O" -U $CT/${F}_single_out_R2.fastq";
fi

if [ -z "$O" ]; then
	echo "No files found for $F"
	exit $E_BADARGS
fi

echo "Bowtie options: $O"

bowtie2 -p 8 -q --no-unal --no-unal -x $HOME/Phage/crAssphage/JQ995537 $O | samtools view -bS - | samtools sort -o $OD/$F.bam
samtools index $OD/$F.bam

