dist<-read.delim("/media/klarqv01/08dcb478-5359-41f4-97c8-469190c8a034/temp/chr1_dist.txt",h=F)
b<-table(cut(dist$V2,breaks=seq(0,2504)))
plot(-1,-1,xlim=c(0,2504),ylim=c(0,ceiling(log10(round(max(b),-6)))))
abline(h=1:7,lwd=2,col="lightgrey")
for(i in 1:6){
  abline(h=i+log10(1:10),col="lightgrey",lty="dashed")
}
points(log10(b+1),pch=20,cex=.5,type="p")
abline(v=which(cumsum(b)>sum(b)*0.8)[1],col="red",lwd=2)
abline(v=which(cumsum(b)>sum(b)*0.9)[1],col="pink",lwd=2)
abline(v=which(cumsum(b)>sum(b)*0.95)[1],col="yellow",lwd=2)
