data.file <- commandArgs(trailingOnly=T)[1]
pdf.file <- commandArgs(trailingOnly=T)[2]
pdf(file=pdf.file, height=3.5, width=8)

par(las=1)
d <- read.table(data.file)

## plot depth
plot.pos <- seq(1,max(d[,2]), by=10)
plot(d[plot.pos,2], d[plot.pos,3], cex=0.1, pch=19, ylab="Depth", xlab="Position on chromosome", main="", xaxt="n")
axis(1,at=c(0,1e6,2e6,3e6,4e6), labels=c("0Mbp","1Mbp", "2Mbp", "3Mbp","4Mbp"), lwd=0.3)

### write text of most common positions
o <- order(d[,3],decreasing=T)
most.common <- d[o,]

text(5e5, c(0.9,0.8,0.7)*most.common[1,3], labels=most.common[1:3,2], cex=0.5, pos=4, offset=0.2)
text(8e5, c(0.9,0.8,0.7)*most.common[1,3], labels=most.common[1:3,3], cex=0.5, pos=4, offset=0.2)
dev.off()
