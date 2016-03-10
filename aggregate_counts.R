# aggregate_counts.R
# Author: Jia Rong Wu
# jwu424 (at) gmail.com 
# 
# Description: Script to aggregate simulated RNA-reads 
# Aggregates each read together by sample and strips the last 5 rows that are only summary statistics on the reads
#
# LICENSE  
# This file is part of the simulate data analysis pipeline.
# aggregate_counts.R is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# any later version.
#
# aggregate_counts.R is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with aggregate_counts.R.  If not, see <http://www.gnu.org/licenses/>.
#

# Open 1 file to get the number of features (genes) prior to loading
x <- read.table("../simulated_reads/htsequencecounts/sample_1_counts.txt", row.names=1, header=F)
n.samples <- nrow(x)	# Number of features (genes)

# Generate the set of sample column names and files to be opened
smpls <- NULL
s.colnames <- NULL
for (i in 1:20)
{
	smpls <- c(smpls,paste("../simulated_reads/htsequencecounts/sample_",i,"_counts.txt", sep=""))
	s.colnames <- c(s.colnames,paste("s_",i,sep=""))
}

# Generate container used to hold the aggregated data
agg.counts <- as.data.frame (matrix(data=NA, nrow=n.samples, ncol = 20))
rownames(agg.counts) <- rownames(x)
colnames(agg.counts) <- s.colnames

# Aggregate Data
for (i in 1:20)
{
	temp <- read.table(smpls[i],row.names=1)
	agg.counts[,i] <- temp	
}

# Strip the last 5 rows from agg.counts
agg.stripped <- agg.counts[1:(n.samples-5),1:ncol(agg.counts)]

# Save data
write.table(agg.counts, file="../aldex_analysis/reads_unadjusted.txt", quote=F, sep="\t")
write.table(agg.stripped, file="../aldex_analysis/reads.txt", quote=F, sep="\t")

# On the agg.stripped... You need to remove the a set and generate 0's into it
# Introduce artificial sparsity to ONE condition
# Last 500 features(genes) of samples in condition B will be "sparse" relative to condition A


# reads_B_500_0.txt
n.samples.adj <- n.samples - 5 
agg.500 <- agg.stripped
agg.500[(n.samples.adj-500):n.samples.adj,11:20] <- 0
write.table(agg.500, file="../aldex_analysis/reads_B_500_0.txt",quote=F,sep="\t")

# reads_A_500_0.txt
agg.500.A.0 <- agg.stripped
agg.500.A.0[(n.samples.adj-500):n.samples.adj,1:10] <- 0
write.table(agg.500.A.0, file="../aldex_analysis/reads_A_500_0.txt",quote=F,sep="\t")



# reads_A_500_MIN.txt
agg.500.A <- agg.500
agg.500.A[(n.samples.adj-500):n.samples.adj,1:10] <- (agg.500.A[(n.samples.adj-500):n.samples.adj,1:10]/100)
agg.500.A <- ceiling(agg.500.A)
write.table(agg.500.A, file="../aldex_analysis/reads_A_500_MIN.txt",quote=F,sep="\t")

# reads_B_500_MIN.txt
agg.500.B.min <- agg.stripped
agg.500.B.min[(n.samples.adj-500):n.samples.adj,11:20] <- (agg.500.B.min[(n.samples.adj-500):n.samples.adj,11:20]/100)
agg.500.B.min <- ceiling(agg.500.B.min)
write.table(agg.500.B.min, file="../aldex_analysis/reads_B_500_MIN.txt",quote=F,sep="\t")



# reads_B_500_ALT.txt
agg.500.B.alt <- agg.stripped
for (i in seq(11,20,2))
{
	agg.500.B.alt[(n.samples.adj-500):n.samples.adj,i] <- agg.500.B.alt[(n.samples.adj-500):n.samples.adj,i] / 1000
}
agg.500.B.alt <- round(agg.500.B.alt, digits = 0)
write.table(agg.500.B.alt, file="../aldex_analysis/reads_B_500_ALT.txt",quote=F,sep="\t")

# reads_A_500_ALT.txt
agg.500.A.alt <- agg.stripped
for (i in seq(1,10,2))
{
	agg.500.A.alt[(n.samples.adj-500):n.samples.adj,i] <- agg.500.A.alt[(n.samples.adj-500):n.samples.adj,i] / 1000
}
agg.500.A.alt <- round(agg.500.A.alt, digits = 0)
write.table(agg.500.A.alt, file="../aldex_analysis/reads_A_500_ALT.txt",quote=F,sep="\t")



# reads_AB_500_ALT.txt
agg.500.AB.alt <- agg.stripped
for (i in seq(1,20,2))
{
	agg.500.AB.alt[(n.samples.adj-500):n.samples.adj,i] <- agg.500.AB.alt[(n.samples.adj-500):n.samples.adj,i] / 1000
}
agg.500.AB.alt <- round(agg.500.AB.alt, digits = 0)
write.table(agg.500.AB.alt, file="../aldex_analysis/reads_AB_500_ALT.txt",quote=F,sep="\t")






# Generate the Meta Table for analysis with other tools 
meta.table <- as.data.frame(matrix(data=NA, nrow=20, ncol=1))
colnames(meta.table) <- "condition"
rownames(meta.table) <- s.colnames 
meta.table[1:10,1] <- 0
meta.table[11:20,1] <- 1
write.table(meta.table, file="../aldex_analysis/meta.txt", quote=F, sep="\t")
