# Specify colour scheme
colors<-paste0(colorRampPalette(c("blue","red"))(10),seq(0,100,length.out = 11))
colors[1]<-paste0(colors[1],"0")
colors[length(colors)]<- substr(colors[length(colors)],1,7)

# Define support functions
plotLDRegion<-function(dataSource, from, to, ...){
  # B is A but sorted for plotting reasons (Z-stack)
  b<-dataSource[dataSource$V3>=from & dataSource$V3 <= to & dataSource$V5 >= from & dataSource$V5 <= to,]
  b<-b[order(b$V12,decreasing = F),]
  plot(b$V3,b$V5,pch=20,cex=.2,col=colors[cut(b$V12,breaks=seq(0,1,length.out = 11),include.lowest = T)],xlim=c(from,to),ylim=c(from,to),xaxs="i",yaxs="i", ...)
}

plotLDRegionTriangular<-function(dataSource, from, to, ...){
  # B is A but sorted for plotting reasons (Z-stack)
  b<-dataSource[dataSource$V3>=from & dataSource$V5<=to & dataSource$V3>=from & dataSource$V5<=to,]
  b<-b[order(b$V12,decreasing = F),]
  plot(b$V3 + ((b$V5-b$V3)/2),b$V5-b$V3,pch=20,cex=.2,col=colors[cut(b$V12,breaks=seq(0,1,length.out = 11),include.lowest = T)],xaxs="i",yaxs="i", ...)
}

# Load some LD data from Tomahawk
ld<-read.delim("1kgp3_chr2_105_1.ld",h=F)
plotLDRegion(ld, 2e6, 5e6, xlab="Coordinates",ylab="Coordinates",main="1KGP3 chr20 2e6-5e6", las=2)
plotLDRegionTriangular(ld, 2e6, 5e6, xlab="Coordinates",ylab="Coordinates",main="1KGP3 chr20 2e6-5e6", las=2)