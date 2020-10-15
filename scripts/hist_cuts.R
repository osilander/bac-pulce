data.file <- commandArgs(trailingOnly=T)[1]
pdf.file <- commandArgs(trailingOnly=T)[2]
pdf(file=pdf.file, height=3.5, width=8)

par(las=1)
du <- read.table(data.file)
d <- du[which(du>0),]

if(length(d)>0) {

  f <- as.data.frame(table(d))
  f[,1] <- as.numeric(as.character(f[,1]))
  ## plot depth
  plot(f[,1], f[,2], ylab="Number of reads", xlab="Position of read start on chromosome", main="", cex=0.2)
  #axis(1,at=c(0,1e6,2e6,3e6,4e6), labels=c("0Mbp","1Mbp", "2Mbp", "3Mbp","4Mbp"), lwd=0.3)

  ### write text of most common positions
  o <- order(f[,2],decreasing=T)
  most.common <- f[o,]

  text(0.08*max(d), c(0.9,0.8,0.7)*most.common[1,2], labels=paste(most.common[1:3,1], "bp", sep=" "), cex=0.5, pos=4, offset=0.2)
  text(0.2*max(d), c(0.9,0.8,0.7)*most.common[1,2], labels=paste(most.common[1:3,2],"reads",sep=" "), cex=0.5, pos=4, offset=0.2)
}
if(length(d)==0) {
  plot(-100,-100,xlim=c(0,1), ylim=c(0,1),xlab="",ylab="",xaxt="n",yaxt="n",main="",bty="n")
 text(0.5,0.5,labels="No mapped reads")
}
dev.off()
