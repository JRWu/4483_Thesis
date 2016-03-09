# Analyze.R
# Author: Jia Rong Wu
# jwu424 (at) gmail.com 
#
# Description: General R script in order to call ALDEx2z (modified) on a simulated dataset
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
# Expected Input Format:
# READS:
#		gene1	gene2	gene...	geneN
#	smpl1
#	smpl2
#	smpl...
#	smplN
#
#
# META:
#		met1	met2	met...	metN
#	smpl1
#	smpl2
#	smpl...
#	smplN
#
#
# Simulated Dataset with fold changes of:
# A		B
# 1 	2*10
# 1 	3*10
# 1		4*10
# 1		5*10
# ...	... (ALL 1 after this point)
# n		1
#
# where "n" refers to the number of features in the whole dataset
# Need to analyze the .called between the two analysis in order to determine how different the two analysis are; this will be the baseline accuracy
# Need to create a set of reduced counts and a set of inflated counts*
# Does the G(x) change between the analysis 1 and 2 between the two transformations?

# Load Required Libraries
library(ALDEx2z)
library(psych)

# Enable Command Line Invocation
args=(commandArgs(trailingOnly=TRUE))
if (length(args)==0)
{
	stop("Expected Usage: Rscript --vanilla Analyze.R filename",call.=FALSE)
} else if (length(args) == 1)
{
	file.name <- as.character(args[1])
}


ald.name <- ("../aldex_analysis/")
rea.file <- (paste(file.name,"_figures/",sep=""))
created.dir <- paste(ald.name,rea.file,sep="")

# Create output directory
dir.create(created.dir)

file.name <- paste(file.name,".txt",sep="")


reads <- read.table(paste("../aldex_analysis/",file.name,sep=""), header=T, row.names=1)
colnames(reads) <- gsub("^X","",colnames(reads))
conds <- c(rep("A",10),rep("B",10))

# Assert each value as numeric 
for (i in 1:ncol(reads))
{
	reads[,i] <- as.numeric(reads[,i])
}

reads <- as.data.frame(reads)

# Aldex with Zero removal on READS
x.clr.z <- aldex.clr(reads,conds,128, TRUE, TRUE, TRUE)
x.e.z <- aldex.effect(x.clr.z, conds, useMC=TRUE) 
x.t.z <- aldex.ttest(x.clr.z, conds)
x.all.z <- data.frame(x.e.z, x.t.z)


# Aldex with NONZERO removal on READS (original)
x.clr.o <- aldex.clr(reads,conds,128, FALSE, TRUE, TRUE)
x.e.o <- aldex.effect(x.clr.o, conds, useMC=TRUE) 
x.t.o  <- aldex.ttest(x.clr.o, conds)
x.all.o <- data.frame(x.e.o, x.t.o)


#HISTOGRAM AND PLOTTING WORK
arrow.len <- 50
# Plot the histograms of default ALDEx2
png(paste(created.dir,"og_btw_cond_hist.png",sep=""),width=800, height=800,res=150)
og.hist <- hist(x.all.o$diff.btw, xlim=c(-10,10), breaks=100, main="Original Btwn-Condition",xlab="Between-Condition Difference", cex.main=2.0, cex.lab=1.5)
og.mean <- mean(x.all.o$diff.btw)
og.med <- median(x.all.o$diff.btw)
og.len <- max(og.hist$counts) 
#og.gmean <- geometric.mean(x.all.o$diff.btw,na.rm=TRUE)
arrows(x0=og.mean, y0=og.len + arrow.len, x1= og.mean, y1 = og.len, length=0.2, col=c(rgb(1,0,0,1)), lwd=4)
arrows(x0=og.med, y0=og.len + arrow.len, x1= og.med, y1 = og.len, length=0.2, col="blue", lwd=4)
#arrows(x0=og.gmean,y0=og.len + arrow.len, x1 = og.gmean, y1 = og.len, length = 0.2, col="green", lwd=4)
abline(v=og.mean, col="red")
abline(v=og.med, col="blue")
#abline(v=og.gmean, col="green")
dev.off()


png(paste(created.dir,"zr_btw_cond_hist.png",sep=""),width=800, height=800,res=150)
# Plot the histogram of Z-adjusted
zr.hist <- hist(x.all.z$diff.btw, xlim=c(-10,10), breaks=100, main="Zero RM Btwn-Condition",xlab="Between-Condition Difference", cex.main=2.0, cex.lab=1.5)
zr.mean <- mean(x.all.z$diff.btw)
zr.med <- median(x.all.z$diff.btw)
zr.len <- max(zr.hist$counts)
#zr.gmean <- geometric.mean(x.all.z$diff.btw)
arrows(x0=zr.mean, y0=zr.len + arrow.len, x1= zr.mean, y1 = zr.len, length=0.2, col=c(rgb(0,0,0,1)), lwd=4)
arrows(x0=zr.med, y0=zr.len + arrow.len, x1= zr.med, y1 = zr.len, length=0.2, col="blue", lwd=4)
#arrows(x0=zr.gmean, y0=zr.len + arrow.len, x1 = zr.gmean, y1 = zr.len, length = 0.2, col="green", lwd=4)
abline(v=zr.mean, col="black")
abline(v=zr.med,col="blue")
#abline(vr=zr.gmean,col="green")
dev.off()

# Set parameters for differences between the two datasets
xlab=NULL
ylab=NULL
xlim=NULL
ylim=NULL
all.col=rgb(0,0,0,0.2)
all.pch=19
all.cex=0.4
called.col="red"
called.pch=20
called.cex=0.6
thres.line.col="darkgrey"
thres.lwd=1.5
test="welch"
cutoff=0.1
rare.col="black"
rare=0
rare.pch=20
rare.cex=0.2


#Determine the "different features"
x <- x.all.o
called <- x$we.eBH <= cutoff
og.called <- as.data.frame(called)
og.names <- which(og.called$called == TRUE)
og.names <- as.data.frame(og.names)
og.ncalled <- length(og.names)

x <- x.all.z
called <- x$we.eBH <= cutoff
zr.called <- as.data.frame(called)
zr.names <- which(zr.called$called == TRUE)
zr.names <- as.data.frame(zr.names)
zr.ncalled <- length(zr.names)

zo.diff <-setdiff(og.names[,1],zr.names[,1])	# Contains the different values 

# Reset to zero removal
x <- x.all.z
png(paste(created.dir,"ALDEx2_1_SIMULATED.png",sep=""),width=800, height=800,res=150)

plot(x$diff.win, x$diff.btw, xlab=xlab, ylab=ylab, col=all.col, pch=all.pch, cex=all.cex, main="ALDEx2.1", cex.main=2.0)
points(x$diff.win[x$rab.all < rare], x$diff.btw[x$rab.all < rare], col=rare.col, pch=rare.pch, cex=rare.cex)
points(x$diff.win[called], x$diff.btw[called], col=called.col, pch=called.pch, cex=called.cex)
abline(0,1, col=thres.line.col, lty=2, lwd=thres.lwd)
abline(0,-1, col=thres.line.col, lty=2, lwd=thres.lwd)
points(x$diff.win[zo.diff], x$diff.btw[zo.diff], col="blue", pch=called.pch, cex=called.cex)
abline(0,0)

dev.off()

x <- x.all.o
called <- x$we.eBH <= cutoff	# Reset called value
png(paste(created.dir,"ALDEx2_0_SIMULATED.png",sep=""),width=800, height=800,res=150)
plot(x$diff.win, x$diff.btw, xlab=xlab, ylab=ylab, col=all.col, pch=all.pch, cex=all.cex, main="ALDEx2.0", cex.main=2.0)
points(x$diff.win[x$rab.all < rare], x$diff.btw[x$rab.all < rare], col=rare.col, pch=rare.pch, cex=rare.cex)
points(x$diff.win[called], x$diff.btw[called], col=called.col, pch=called.pch, cex=called.cex)
abline(0,1, col=thres.line.col, lty=2, lwd=thres.lwd)
abline(0,-1, col=thres.line.col, lty=2, lwd=thres.lwd)
abline(0,0)
dev.off()


# Difference Between, Difference Within, Effects
pdf(paste(created.dir,"Zero_CLR_vs_Org_CLR_Diff_BTW.pdf",sep=""))
plot(x.e.o$diff.btw, x.e.z$diff.btw, xlab="Original CLR Diff Btw", ylab="Zero CLR Diff.Btw", main="Zero Diff Between vs CLR Diff Between", pch=19, col=rgb(0,0,0,0.2))
points(x.e.o$diff.btw[zo.diff], x.e.z$diff.btw[zo.diff],col="blue",pch=20)
abline(0,1)
dev.off()

pdf(paste(created.dir,"Zero_CLR_vs_Org_CLR_Diff_WIN.pdf",sep=""))
plot(x.e.o$diff.win, x.e.z$diff.win, xlab="Original CLR Diff Win", ylab="Zero CLR Diff.Win", main="Zero Diff Within vs CLR Diff Within", pch=19, col=rgb(0,0,0,0.2))
points(x.e.o$diff.win[zo.diff], x.e.z$diff.win[zo.diff],col="blue",pch=20)
abline(0,1)
dev.off()

pdf(paste(created.dir,"Zero_CLR_vs_Org_CLR_Effects.pdf",sep=""))
plot(x.e.o$effect, x.e.z$effect, xlab="Original CLR Effect", ylab="Zero CLR Effect", main="Zero Effect vs CLR Effect", pch=19, col=rgb(0,0,0,0.2))
points(x.e.o$effect[zo.diff], x.e.z$effect[zo.diff],col="blue",pch=20)
abline(0,1)
dev.off()

# Original aldex plotting Graphs
pdf(paste(created.dir,"CLR_ORIGINAL.pdf",sep=""))
aldex.plot(x.all.o)
abline(0,0)
dev.off()

pdf(paste(created.dir,"CLR_ZERO_REMOVAL.pdf",sep=""))
aldex.plot(x.all.z)
abline(0,0)
dev.off()

#Features that are different between the transformations
n.diff <- length(zo.diff)
print(n.diff)

# Summary stats MIN
og.diff.btw.min <- min(x.e.o$diff.btw[zo.diff])
zr.diff.btw.min <- min(x.e.z$diff.btw[zo.diff])
og.diff.win.min <- min(x.e.o$diff.win[zo.diff])
zr.diff.win.min <- min(x.e.z$diff.win[zo.diff])
og.effect.min <- min(x.e.o$effect[zo.diff])
zr.effect.min <- min(x.e.z$effect[zo.diff])

# Summary stats MAX
og.diff.btw.max <- max(x.e.o$diff.btw[zo.diff])
zr.diff.btw.max <- max(x.e.z$diff.btw[zo.diff])
og.diff.win.max <- max(x.e.o$diff.win[zo.diff])
zr.diff.win.max <- max(x.e.z$diff.win[zo.diff])
og.effect.max <- max(x.e.o$effect[zo.diff])
zr.effect.max <- max(x.e.z$effect[zo.diff])

# Summary stats MEAN
og.diff.btw.mean <- mean(x.e.o$diff.btw[zo.diff])
zr.diff.btw.mean <- mean(x.e.z$diff.btw[zo.diff])
og.diff.win.mean <- mean(x.e.o$diff.win[zo.diff])
zr.diff.win.mean <- mean(x.e.z$diff.win[zo.diff])
og.effect.mean <- mean(x.e.o$effect[zo.diff])
zr.effect.mean <- mean(x.e.z$effect[zo.diff])

# Summary stats MEDIAN
og.diff.btw.median <- median(x.e.o$diff.btw[zo.diff])
zr.diff.btw.median <- median(x.e.z$diff.btw[zo.diff])
og.diff.win.median <- median(x.e.o$diff.win[zo.diff])
zr.diff.win.median <- median(x.e.z$diff.win[zo.diff])
og.effect.median <- median(x.e.o$effect[zo.diff])
zr.effect.median <- median(x.e.z$effect[zo.diff])


#Generate a table of summary statistics
summ.rownames <- c("D.Btw","D.Win","Effect")
summ.colnames <- c("og.Min", "zr.Min", "og.Max", "zr.Min", "og.Mean", "zr.Mean", "og.Median", "zr.Median","og.nCalled","zr.nCalled")
summ.stats <- matrix(data=NA, nrow=3, ncol=10)
rownames(summ.stats) <- summ.rownames
colnames(summ.stats) <- summ.colnames

summ.stats[1,1] <-og.diff.btw.min
summ.stats[1,2] <-zr.diff.btw.min
summ.stats[1,3] <-og.diff.btw.max
summ.stats[1,4] <-zr.diff.btw.max
summ.stats[1,5] <-og.diff.btw.mean
summ.stats[1,6] <-zr.diff.btw.mean
summ.stats[1,7] <-og.diff.btw.median
summ.stats[1,8] <-zr.diff.btw.median
summ.stats[1,9] <- og.ncalled
summ.stats[1,10] <- zr.ncalled

summ.stats[2,1] <-og.diff.win.min
summ.stats[2,2] <-zr.diff.win.min
summ.stats[2,3] <-og.diff.win.max
summ.stats[2,4] <-zr.diff.win.max
summ.stats[2,5] <-og.diff.win.mean
summ.stats[2,6] <-zr.diff.win.mean
summ.stats[2,7] <-og.diff.win.median
summ.stats[2,8] <-zr.diff.win.median
summ.stats[2,9] <- ""
summ.stats[2,10] <- ""

summ.stats[3,1] <-og.effect.min
summ.stats[3,2] <-zr.effect.min
summ.stats[3,3] <-og.effect.max
summ.stats[3,4] <-zr.effect.max
summ.stats[3,5] <-og.effect.mean
summ.stats[3,6] <-zr.effect.mean
summ.stats[3,7] <-og.effect.median
summ.stats[3,8] <-zr.effect.median
summ.stats[3,9] <- ""
summ.stats[3,10] <- ""

write.table(summ.stats, file=paste(created.dir,n.diff,"_summary_stats_SIMULATED.txt",sep=""), quote=F, sep="\t")
