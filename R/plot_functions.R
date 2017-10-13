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
plotD<-function(dataSource, min_snps = 5, ...){
  # Store plot parameters
  temp<-par()
  
  # Filter data
  b<-dataSource[dataSource$n_snps>=min_snps,]
  
  # Plot  
  par(mfrow=c(5,1))
  par(mar=c(0, 5, 3, 3))
  plot(b$cumBinFrom,1/(b$n_snps/(b$cumBinTo-b$cumBinFrom)),pch=20,las=2,xaxt='n',xaxs='i',ylab="SNV/bp",...)
  abline(v=d[which(!duplicated(b$contigID)),"cumBinFrom"],lty="dashed",col="grey")
  par(mar=c(0, 5, 0, 3))
  plot(b$cumBinFrom,-log10(b$meanMAF),pch=20,las=2,xaxt='n',xaxs='i',ylab="-log10(mean MAF)",...)
  abline(v=d[which(!duplicated(b$contigID)),"cumBinFrom"],lty="dashed",col="grey")
  
  par(mar=c(0, 5, 0, 3))
  plot(b$cumBinFrom,-log10(b$het),pch=20,las=2,xaxt='n',xaxs='i',ylab="-log10(mean het)",...)
  abline(v=d[which(!duplicated(b$contigID)),"cumBinFrom"],lty="dashed",col="grey")
  
  par(mar=c(0, 5, 0, 3))
  plot(-1,-1,xlim=c(min(b$cumBinFrom),max(b$cumBinFrom)), ylim=c(min(b$TajimaD),max(b$TajimaD)), las=2,xaxt='n',xaxs='i',ylab="Tajima's D",...)
  rect(min(b$cumBinFrom),-2,max(b$cumBinFrom), 2, col="#F5DEB350",border = NA);
  abline(h=0,lwd=3,col="lightgrey")
  points(b$cumBinFrom,b$TajimaD,pch=20, ...)
  abline(v=d[which(!duplicated(b$contigID)),"cumBinFrom"],lty="dashed",col="grey")
  par(mar=c(3, 5, 0, 3))
  plot(b$cumBinFrom,b$meanPI,pch=20,las=2,ylab="mean Pi",xaxs='i',xlab="Cumulative genomic position",...)
  abline(v=d[which(!duplicated(b$contigID)),"cumBinFrom"],lty="dashed",col="grey")
  # Restore plot parameters
  par(mfrow=temp$mfrow,mar=temp$mar)
}

d<-read.delim("~/Desktop/1000GP3/tajima_1kgp3.txt",header = T)
plotD(d[d$n_snps>50,],cex=.1)


# temp
plot(-1,-1,ylim=c(min(d2[,-c(1:3,ncol(d2))]),max(d2[,-c(1:3,ncol(d2))])),xlim=c(min(d2$V2),max(d2$V2)))
