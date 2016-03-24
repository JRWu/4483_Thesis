# simulate_data.R
# Author: Jia Rong Wu
# jwu424 (at) gmail.com 
# 
# Description: Script that utilizes the polyester package in order to simulate RNA-seq reads
# 
# LICENSE  
# This file is part of the simulate data analysis pipeline.
# Analyze.R is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# any later version.
#
# Analyze.R is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with aggregate_counts.R.  If not, see <http://www.gnu.org/licenses/>.
#
# uses complete_genome.fna from yeast
# seq/barton_RNAseq/yeast
# /bowtie2-2.2.7/bowtie2-build ../source_files_barton_yeast/complete_genome.fna index/yeast_index
# mapping_script.sh will take these reads and then map them to the yeast genome
#
# For each file in the yeast directory 
# .asn		Appears to be "metadata"
# .faa		Appears to be the "amino acid sequence" of gene product
# .ffn	*	Appears to be complete sequence of gene
# .fna		(Fasta Nucleic Acid) Appears to be the "complete sequence" 
# .frn		Appears to be predicted sequences
# .gbk		More Metadata
# .gff	*	Chromosome information (yeast.gff contains each chromosome)
# .ptt	 	More Sequence Region
# .rnt		Even More Sequence region
# .rpt		Metadata
# .val		Cannot read format
#
# To get the # of lines, use
# wc -l in the yeast.gff file
# gives you 6349 lines
# Need to strip 6349 files from the set of values
# Sample_n_counts.txt contains tab-delimited values
# 

library(polyester)
library(Biostrings)

fasta = readDNAStringSet("source_files_barton_yeast/yeast_grouped.fa")

full_fasta = fasta

# Reads per transcript = transcriptLength / readLength * coverage
# FPKM: fragments per kilobase of exon per million reads mapped
coverage <- 20	# Num of times a genome has been sequenced (depth)
readspertx = round(coverage * width(full_fasta) / 100)


# Ask greg what a realistic fold change expectation would be in a dataset
# What would be a good fold change value for the "ideal" dataset
# Generate a list of fold changes, There is 10x range(2,4) fold change between the two

# So you want to discuss the fold change between the two conditions
# Use replicates that are CLOSE to each other, and some replicates that are NOT close to each other 
# Since the effect size is what's being examined in ALDEx, effect = btw/win
# Win will be driving the difference 

fold_changes = matrix(c(rep(1,length(fasta)),rep(2,10),rep(3,10),rep(4,10),rep(5,10), rep(1,(length(fasta)) - 40)       ), nrow = length(fasta))

# simulate_exp writes files to the driectory the file is in 
# num_reps = 10 + 10
# Samples 1-10 = condition 1 (Group A)
# Samples 11-20 = condition 2 (Group B)
simulate_experiment("source_files_barton_yeast/ERCC92.fa", reads_per_transcript=readspertx, num_reps=c(10,10), fold_changes=fold_changes, outdir='simulated_reads', paired=FALSE) 

# Need to call the mapping_script to assemble it together
# /bowtie2-2.2.7 -p 8 -x 



dir.create("simulated_reads/fastq_files")
dir.create("simulated_reads/samfiles")
dir.create("simulated_reads/mapped")
dir.create("simulated_reads/mapped_fixed")
dir.create("simulated_reads/htsequencecounts")
