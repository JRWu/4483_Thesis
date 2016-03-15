# PRIOR TO EXECUTION
# Extracted yeast ffn files and appended them together 
# Generated a bowtie index the yeast genome
# Mapped the aggregate of ffn files to the complete yeast genome with bowtie through mapping.sh

# EXECUTION
# Handle the calling of the simulate_data.R
# 	This generates the directories required for analysis afterwards
# 	Handle making directories called fastq_files and samfiles
# Run fasta_to_fq in order to add dummy quality scores for the alignment
# Runs mapping.sh to create the samfiles from the set of fq files
# Invokes mapping_script.sh to determine the counts from the SAMfiles
# Calls aggregate_counts.R in order to append all the data together from the separate sample countfiles
# Calls Analyze.R for ALDEx2.1 analysis on the dataset 
	# This allows for generation of graphs and etc to determine differences  
	# Displays the histograms and various data 
	
#R CMD BATCH ../simulate_data.R
Rscript --vanilla Analyze.R reads # replace reads.txt with your filename
Rscript --vanilla Analyze.R reads_A_500_0
Rscript --vanilla Analyze.R reads_B_500_0

Rscript --vanilla Analyze.R reads_A_500_MIN
Rscript --vanilla Analyze.R reads_B_500_MIN

Rscript --vanilla Analyze.R reads_A_500_ALT
Rscript --vanilla Analyze.R reads_B_500_ALT

Rscript --vanilla Analyze.R reads_AB_500_ALT

./fasta_to_fq.sh 
./mapping.sh
./mapping_script.sh
R CMD BATCH aggregate_counts.R
R CMD BATCH Analyze.R 

# Some magic with saving the geometric means happens here

Rscript --vanilla plot_strips.R