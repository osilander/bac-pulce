data.file <- commandArgs(trailingOnly=T)[1]
pdf.file <- commandArgs(trailingOnly=T)[2]
pdf(file=pdf.file, height=3.5, width=3.5)

par(las=1)
du <- read.table(data.file)
d <- du[which(du>0),1]

d <- du[which(du>0),1]

if(length(d)>0) {
  ## plot cumulative
  plot(d, 1:length(d)/length(d), ty="l", ylab="Cumulative fraction of reads", xlab="Position of read start on chromosome", main="", cex=0.2)
  abline(h=seq(0,1,by=0.05), lty=2, lwd=0.2) 
}

if(length(d)==0) {
  plot(-100,-100,xlim=c(0,1), ylim=c(0,1),xlab="",ylab="",xaxt="n",yaxt="n",main="",bty="n")
 text(0.5,0.5,labels="No mapped reads")
}

dev.off()
