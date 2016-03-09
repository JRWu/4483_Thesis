# mapping.sh
# Author: Gloor-lab Western University
# Modified by: Jia Rong Wu
# jwu424 (at) gmail.com
#
# Description: Script that maps the respective fq files to the original genome to generate samfiles
# Build index
# $BIN/bowtie2-build complete_genome.fna index/complete

INDX="index/yeast_index"
FAQLOCATION="../simulated_reads/fastq_files/*"

SIMRDS="../simulated_reads"
FQFILES="/fastq_files/"
FQEXT=".fq"

SMFILES="/samfiles/"
SM=".sam"

for file in $FAQLOCATION
do
	# Get the name of the file only without the folders, but still has extension
	filename=$(basename $file)
	# Get rid of the .fq
	filename="${filename%.*}"
	# Execute the command that's the whole point of this
	# Edited "$BIN/bowtie2" to "$BIN"bowtie2"" because the original was writing ...//bowtie2. Not sure if that matters
	echo $filename	
	/bowtie2-2.2.7/bowtie2 -x $INDX -U $SIMRDS$FQFILES$filename$FQEXT -p 8 -S $SIMRDS$SMFILES$filename$SM
done
