# 
# Calculate the number of positions within a region with no coverage.


if [ -z $1 ]; then
	echo "BAM file to check?"
	exit
fi

A=$(samtools depth -a -r "JQ995537:25633-26964" $1  | grep -e '\s0$' | wc -l);
B=$(samtools depth -a -r "JQ995537:33708-35062" $1  | grep -e '\s0$' | wc -l);
C=$(samtools depth -a -r "JQ995537:43819-45057" $1  | grep -e '\s0$' | wc -l);


echo -e "$1\t$A\t$B\t$C";
