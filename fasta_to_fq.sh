COUNTER=0

SMP="sample_"
FASTA=".fasta"
UN=".fq"
DX="2.fq"
UP="../"
SIMR="simulated_reads/"
FAQF="fastq_files/"
OH="0"

# 2 separate loops because the polyester package appends 0 in front of names in the 10's digits
# 1-9
while [ $COUNTER -lt 9 ]; do
	echo $COUNTER
	let COUNTER=COUNTER+1
	perl fasta_to_fastq.pl $UP$SIMR$SMP$OH$COUNTER$FASTA > $UP$SIMR$FAQF$SMP$COUNTER$UN
	echo $UP$SIMR$FAQF$SMP$COUNTER$UN
done

#10-"n"
while [ $COUNTER -lt 20 ]; do
	echo $COUNTER
	let COUNTER=COUNTER+1
	perl fasta_to_fastq.pl $UP$SIMR$SMP$COUNTER$FASTA > $UP$SIMR$FAQF$SMP$COUNTER$UN
	echo $UP$SIMR$FAQF$SMP$COUNTER$UN
done
