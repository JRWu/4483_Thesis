x <- read.table("sample_gmeans_orig.txt", header=T, row.names=1)
y <- read.table("sample_gmeans_zero.txt", header=T, row.names=1)

maxim <- max(x,y)
minim <- min(x,y)

diff <- abs(maxim-minim)

zgmeans <- y
ogmeans <- x

zgmeans.mean <- t(as.data.frame(colMeans(zgmeans)))
ogmeans.mean <- t(as.data.frame(colMeans(ogmeans)))

z.a <- zgmeans.mean[1:10]
z.b <- zgmeans.mean[11:20]
o.a <- ogmeans.mean[1:10]
o.b <- ogmeans.mean[11:20]

zr.diff <- 2^(abs(mean(z.a)-mean(z.b)))
og.diff <- 2^(abs(mean(o.b)-mean(o.a)))
