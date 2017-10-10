colors<-paste0(colorRampPalette(c("blue","red"))(10),seq(0,100,length.out = 11))
colors[1]<-paste0(colors[1],"0")
colors[length(colors)]<- substr(colors[length(colors)],1,7)

plotLDRegion<-function(dataSource, ...){
  b<-dataSource[order(dataSource$R2,decreasing = F), ]
  plot(b$Aposition,b$Bposition,
       pch=20,cex=.2,
       col=colors[cut(b$R2,breaks=seq(0,1,length.out = 11),include.lowest = T)],
       xaxs="i",yaxs="i", 
       ...
  )
}

# Tajima D
d<-read.delim("~/Desktop/cichlids__tajima.txt",skip = 1,header = F)
par(mfrow=c(3,1))
plot(d$V5,-log10(d$V9),pch=20,cex=.3)
plot(d$V5,d$V8,pch=20,cex=.3)
plot(d$V5,d$V7,pch=20,cex=.3)

plotD<-function(dataSource, ...){
  par(mfrow=c(4,1))
  par(mar=c(0, 5, 3, 3))
  plot(dataSource$cumBinFrom,-log10(dataSource$meanMAF),pch=20,las=2,xaxt='n',ylab="mean MAF",...)
  par(mar=c(0, 5, 0, 3))
  plot(dataSource$cumBinFrom,dataSource$TajimaD,pch=20,las=2,xaxt='n',ylab="Tajima's D",...)
  abline(h=0,lwd=3,col="lightgrey")
  par(mar=c(0, 5, 0, 3))
  plot(dataSource$cumBinFrom,dataSource$meanPI,pch=20,las=2,xaxt='n',ylab="mean Pi",...)
  par(mar=c(3, 5, 0, 3))
  plot(dataSource$cumBinFrom,dataSource$k_hat,pch=20,las=2,ylab="k_hat",xlab="Cumulative genomic position",...)
}