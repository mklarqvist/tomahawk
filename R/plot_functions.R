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