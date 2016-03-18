# Analyze.R
# Author: Jia Rong Wu
# jwu424 (at) gmail.com 
#
# Description: R script to plot the geometric means of each sample based on ALDEx2 or ALDEx2z analysis
#
# LICENSE  
# This file is part of the simulate data analysis pipeline.
# Analyze.R is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# any later version.
#
# plot_strips.R is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with plot_strips.R.  If not, see <http://www.gnu.org/licenses/>.

# No External packages required
# Set the labels for reading in 
ald.name <- ("../aldex_analysis/")
nom <- c("reads_A_500_0","reads_A_500_ALT","reads_A_500_MIN","reads_AB_500_ALT","reads_B_500_0","reads_B_500_ALT","reads_B_500_MIN","reads")
fg <- "_figures/"
znom <- "sample_gmeans_zero.txt"
onom <- "sample_gmeans_orig.txt"

mpdf <- "stripchart_mean.pdf"
ipdf <- "stripchart_instance.pdf"

dname <- "GMean_density.png"



# Iterate and plot stripchart of each dataset
for (i in 1:8)
{
	finz <- paste(ald.name,nom[i],fg,znom,sep="")
	fino <- paste(ald.name,nom[i],fg,onom,sep="")

	pdfm <- paste(ald.name,nom[i],fg,mpdf,sep="")
	pdfi <- paste(ald.name,nom[i],fg,ipdf,sep="")

	pngz <- paste(ald.name,nom[i],fg,dname,sep="")

	# Set the name of the analysis
	manlbl=paste(nom[i]," ", sep="")

	#Read Tables		
	zgmeans <- read.table(finz, header=T, row.names=1, sep="\t")	
	ogmeans <- read.table(fino, header=T, row.names=1, sep="\t")

	# Transpose tables and get the average geometric mean per sample
	zgmeans.mean <- t(as.data.frame(colMeans(zgmeans)))
	ogmeans.mean <- t(as.data.frame(colMeans(ogmeans)))
	
	# Empty matrix to store the set of geometric means
	z.new <- matrix(, nrow=20, ncol=16)
	o.new <- matrix(, nrow=20, ncol=16)
	
	z.new[,1] <- zgmeans.mean	# Save Means
	o.new[,1] <- ogmeans.mean	# Save Means
	
	# Append dummy data because means are already saved
	# Stripcharts won't plot 20 samples if there aren't 20 rows
	for (i in 2: 16)
	{
		z.new[,i] <- 0
		o.new[,i] <- 0
	}
	z.new <- as.data.frame(t(z.new))
	o.new <- as.data.frame(t(o.new))	

	# Save the names of the datasets	
	colnames(z.new) <- colnames(zgmeans)
	rownames(z.new) <- rownames(zgmeans)
	colnames(o.new) <- colnames(ogmeans)
	rownames(o.new) <- rownames(ogmeans)

	# Specific range for plots
	ymax <- max(zgmeans,ogmeans)	# Positive
	ymin <- min(zgmeans,ogmeans)	# Negative

	pdf(pdfm)
	stripchart(z.new,xlim=c(ymax,ymin), col=c(rgb(0,0,1,0.7)), pch=15, xlab=("Geometric Mean G(x)"), ylab=("Sample"),main=paste(manlbl,"Mean",sep=""))	# Blue is zero removed
	stripchart(o.new,add=TRUE,col=c(rgb(1,0,0,0.7)),pch=15)	# Red is original
	dev.off()

	pdf(pdfi)
	stripchart(zgmeans,method="jitter",jitter=2,xlim=c(ymax,ymin), col=c(rgb(0,0,1,0.4)), pch=15, xlab=("Geometric Mean G(x)"), ylab=("Sample"),main=paste(manlbl,"Instance", sep=""))	# Blue is zero removed
stripchart(ogmeans,method="jitter",jitter=2,add=TRUE,col=c(rgb(1,0,0,0.4)),pch=15)	# Red is original
	dev.off()
	
	density.A.z <- NULL
	density.B.z <- NULL
	density.A.o <- NULL
	density.B.o <- NULL
	
	for (k in 1:10)
	{
		density.A.z <- c(density.A.z, zgmeans[,k])
		density.A.o <- c(density.A.o, ogmeans[,k])
	}
	for (k in 11:20)
	{
		density.B.z <- c(density.B.z, zgmeans[,k])
		density.B.o <- c(density.B.o, ogmeans[,k])
	}
	density.A.z <- density(density.A.z)
	density.B.z <- density(density.B.z)
	
	density.A.o <- density(density.A.o)
	density.B.o <- density(density.B.o)
	
	xmx <- max(density.A.z$y, density.B.z$y, density.A.o$y, density.B.o$y)
	xmn <- 0
	thickn <- 2
	png(pngz, width=4,height=4,units="in",res=150)
	plot (density.A.z,xlim=c(ymax,ymin), ylim=c(xmn,xmx), col="Blue", main=paste(manlbl, "Density", sep=" "), lwd=thickn, xlab="Geometric Mean G(x)", xaxp=c(round(ymax,2),round(ymin,2),10))
	lines(density.B.z, col="Blue", lwd=thickn)
	lines(density.A.o, col="Red", lwd=thickn)
	lines(density.B.o, col="Red", lwd=thickn)
	dev.off()
	
}



# Generate the binary heatmap figures
tx <- ".txt"

for (i in 1:8)
{
	heatname <- paste(ald.name,nom[i],tx,sep="")
	heatout <- paste(ald.name,nom[i],fg,"heatmap.pdf",sep="")
	pdf(heatout)
	x <- read.table(heatname, header=T, row.names=1)
	rownames(x) <- NULL
	x[x>1] <- 1
	heatmap(as.matrix(x[5700:6000,]), Rowv=NA,Colv=NA, col=c("white","black"),scale="column",labRow="")
	dev.off()
}




