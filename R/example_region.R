# Specify colour scheme
colors<-paste0(colorRampPalette(c("blue","red"))(10),seq(0,100,length.out = 11))
colors[1]<-paste0(colors[1],"0")
colors[length(colors)]<- substr(colors[length(colors)],1,7)

# Define support functions
plotLDRegion<-function(dataSource, from, to, upper = FALSE, lower = FALSE, add = FALSE, ...){
  # Assumes all the data is from the same chromosome
  if(upper == TRUE){
    b<-dataSource[dataSource$POS_A >= from & dataSource$POS_A <= to & dataSource$POS_B >= from & dataSource$POS_B <= to & dataSource$POS_A < dataSource$POS_B,]
  } else if(lower == TRUE){
    b<-dataSource[dataSource$POS_A >= from & dataSource$POS_A <= to & dataSource$POS_B >= from & dataSource$POS_B <= to & dataSource$POS_B < dataSource$POS_A,]
  } else {
    b<-dataSource[dataSource$POS_A >= from & dataSource$POS_A <= to & dataSource$POS_B >= from & dataSource$POS_B <= to,]
  }
  b$R2[b$R2>1]<-1 # In cases of rounding error
  b<-b[order(b$R2,decreasing = F),] # sort for Z-stack
  if(add == TRUE){
    points(b$POS_A,b$POS_B,pch=20,cex=.2,col=colors[cut(b$R2,breaks=seq(0,1,length.out = 11),include.lowest = T)], ...)
  } else {
    plot(b$POS_A,b$POS_B,pch=20,cex=.2,col=colors[cut(b$R2,breaks=seq(0,1,length.out = 11),include.lowest = T)],xlim=c(from,to),ylim=c(from,to),xaxs="i",yaxs="i", ...)
    abline(0,1,lwd=2,col="grey")
  }
}

plotLDRegionTriangular<-function(dataSource, from, to, ...){
  # Assumes all the data is from the same chromosome
  b<-dataSource[dataSource$POS_A>=from & dataSource$POS_A<=to & dataSource$POS_B>=from & dataSource$POS_B<=to,]
  b<-b[b$POS_A<b$POS_B,] # upper triangular only
  b<-b[order(b$R2,decreasing = F),] # sort for Z-stack
  plot(b$POS_A + ((b$POS_B-b$POS_A)/2),b$POS_B-b$POS_A,pch=20,cex=.2,col=colors[cut(b$R2,breaks=seq(0,1,length.out = 11),include.lowest = T)],xaxs="i",yaxs="i", ...)
}

# Load some LD data from Tomahawk
ld<-read.delim("1kgp3_chr2_105_1.ld",h=F)
plotLDRegion(ld, 2e6, 5e6, xlab="Coordinates",ylab="Coordinates",main="1KGP3 chr20 2e6-5e6", las=2)
plotLDRegionTriangular(ld, 2e6, 5e6, xlab="Coordinates",ylab="Coordinates",main="1KGP3 chr20 2e6-5e6", las=2)