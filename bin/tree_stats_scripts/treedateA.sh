export PRIMER=A
perl ../../bin/treestats.pl -t Primer$PRIMER/seqs.$PRIMER.rc.trim.aln.treefile -m Primer$PRIMER/metadata.$PRIMER.txt -n date -r 5 > treestats$PRIMER/date.$PRIMER.$SGE_TASK_ID.txt
