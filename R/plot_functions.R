colors<-paste0(colorRampPalette(c("blue","red"))(10),seq(0,100,length.out = 11))
colors[1]<-paste0(colors[1],"0")
colors[length(colors)]<- substr(colors[length(colors)],1,7)

occurences<-sort(table(ld$V4),decreasing = T)
multiples<-as.numeric(names(occurences[occurences>5]))
#Decay
decay<-function(pos, ...){
  par(mar=c(2,2,2,2))
  plot(ld[ld$V4==multiples[pos],6] - ld[ld$V4==multiples[pos],4],ld[ld$V4==multiples[pos],13],pch=20,ylim=c(0,1),cex=(-10*(log10(.5)+log10(.5)))/ld[ld$V4==multiples[pos],2], ...)
}
decay(1)

##test
testModel<-function(pos){
  dat<-data.frame("y"=a[a$V3==multiples[pos],12], 
                  "x"=a[a$V3==multiples[pos],5])
  dat$x<-dat$x-min(dat$x)+1
  
  #plot(dat$x,dat$y)
  mod <- nls(y ~ a*x^(-a*b), data = dat, start = list(a = 1, b = 0.15),algorithm="port",weights = a[a$V3==multiples[pos],12])
  # plot decay
  modelRatio<-(coef(mod)["a"]*dat$x^(-coef(mod)["a"]*coef(mod)["b"]))/mean(dat$y)
  
  par(mfrow=c(2,1))
  decay(pos,col=c("blue","red")[as.factor(modelRatio>2)])
  lines(a[a$V3==multiples[pos],5],coef(mod)["a"]*1/dat$x^coef(mod)["b"],lwd=2,col="red")
  abline(h=mean(dat$y),lwd=2,col="blue",lty="dashed")
  plot(modelRatio,type="l")
}

plot(a$V3,a$V5,pch=20,cex=.5,col=colors[cut(a$V12,breaks=seq(0,1,length.out = 11),include.lowest = T)],xlim=c(750e3,950e3),ylim=c(750e3,950e3))

#
estimator<-function(pos){
  a1<- 1 - pnorm((a[a$V3==multiples[pos],12]-mean(a[a$V3==multiples[pos],12]))/sd(a[a$V3==multiples[pos],12]),lower.tail = F)
  b1<- -log10(a[a$V3==multiples[pos],13])
  #plot(a1)
  composite<- a1 * a[a$V3==multiples[pos],12]
  plot(composite,ylim=c(0,1))
  return(composite)
}
x<-seq(mean(a[a$V3==multiples[1],12])-3*sd(a[a$V3==multiples[1],12]),mean(a[a$V3==multiples[1],12])+3*sd(a[a$V3==multiples[1],12]),length=1000)
plot(x,dnorm(x,mean(a[a$V3==multiples[1],12]),sd(a[a$V3==multiples[1],12])),type="l")

gtf<-read.delim("~/Documents/Homo_sapiens.GRCh37.75_main.txt",head=F)

from<-min(ld$V4[ld$V3==0&ld$V5==0])
to<-max(ld$V6[ld$V3==0&ld$V5==0])

plotLDRegion<-function(dataSource, from, to, ...){
  #temp<-gtf[(gtf$V4<=from&gtf$V5>from)|(gtf$V4>from&gtf$V4<to)|(gtf$V4>=from&gtf$V5<=to)|(gtf$V5>from&gtf$V5<=to),]
  #temp<-temp[temp$V2=="protein_coding",]
  #layout(mat = c(1,2,3),heights = c(1,2,8))
  #par(mar=c(0,3,3,3))
  #plot(-1,-1,ylim=c(0,1),xlim=c(min(temp$V4),max(temp$V5)),xaxt="n",xaxs="i")
  #rect(temp$V4,0.1,temp$V5,0.9,col="black")
  #par(mar=c(0,3,0,3))
  #plot(seq(from,to,by=1000)[-1],table(cut(dataSource[dataSource$V4>=from&dataSource$V6<=to,3],seq(from,to,by=1000))),pch=20,cex=1,xaxt="n",xaxs="i")
  #par(mar=c(3,3,0,3))
  # B is A but sorted for plotting reasons (Z-stack)
  b<-dataSource[dataSource$V2<30,]
  b<-b[order(b$V13,decreasing = F),]
  plot(b$V4,b$V6,pch=20,cex=.2,col=colors[cut(b$V13,breaks=seq(0,1,length.out = 11),include.lowest = T)],xlim=c(from,to),ylim=c(from,to),xaxs="i",yaxs="i", ...)
}

for(i in 1:49){ 
  from=(i-1)*5e6; 
  to=(i*5e6);
  filename = sprintf("~/Desktop/chr1_slices/metabric/chr1_block%i.jpeg",i)
  jpeg(filename = filename,width = 2250,height = 1500, pointsize = 10,units = "px")
  #temp<-gtf[gtf$V1==1&((gtf$V4<=from&gtf$V5>from)|(gtf$V4>from&gtf$V4<to)|(gtf$V4>=from&gtf$V5<=to)|(gtf$V5>from&gtf$V5<=to)),]
  #temp<-temp[temp$V2=="protein_coding"&temp$V3=="gene",]
  #layout(mat = c(1,2,3),heights = c(1,2,8))
  #par(mar=c(0,3,3,3))
  #plot(-1,-1,ylim=c(0,1),xlim=c(from, to),xaxt="n",xaxs="i")
  #if(nrow(temp) != 0)
  #  rect(temp$V4,0.1,temp$V5,0.9,col="black")
  #par(mar=c(0,3,0,3))
  #plot(seq(from,to,by=1000)[-1],table(cut(ld[ld$V4>=from&ld$V4<=to&ld$V6>=from&ld$V6<=to,4],seq(from,to,by=1000))),pch=20,cex=1,xaxt="n",xaxs="i")
  
  #layout(mat = c(1,2),heights = c(1,4))
  #par(mar=c(0,3,1,3))
  plot.new()
  pushViewport(viewport(y = 0.6, height = 0.2, width = 0.8, angle=135))
  plotRecomb(recomb, "chr1", from, to, xaxt='n', las=2, newpage=FALSE)
  #par(mar=c(0,3,0,3))
  #segRegion<-seg[seg$Chromosome==1&((seg$Start>=from&seg$Start<=to)|(seg$End>=from&seg$End<=to)),]
  #segRegion<-segRegion[abs(segRegion$Segment_Mean)>0.5,]
  recombTemp<-recomb[recomb$V3>=from&recomb$V4<=to&recomb$V2=="chr1",]
  #plot(-1,-1,xlim=c(from,to),ylim=c(min(segRegion$Segment_Mean),max(segRegion$Segment_Mean)),xaxs="i")
  #rect(from, -0.5, to, 0.5, col="lightgrey",border=NA)
  #points(segRegion$Start,segRegion$Segment_Mean,pch=20)
  #points(segRegion$End,segRegion$Segment_Mean,pch=17)
  #abline(v=recomb$V3[recomb$V2=="chr1"],col="grey",lty="dashed")
  par(mar=c(3,3,0,3))
  internal <- ld[ld$V4>=from-(to-from)&ld$V4<=to&ld$V6>=from&ld$V6<=to+(to-from),c(4,6,13)]
  plot(internal$V4+internal$V6,internal$V6-internal$V4,pch=20,cex=.25,col=colors[cut(internal$V13,breaks=seq(0,1,length.out = 11),include.lowest = T)],xlim=c(from*2,to*2),ylim=c(0, 500e3),xaxs="i",yaxs="i",las=2)
  abline(v=recomb$V3[recomb$V2=="chr1"],col="grey",lty="dashed")
  dev.off()
}

plot(unlist(lapply(split(test,test$V3),function(a)sum(a$V12>0.5))),pch=20,type="l")

scaleFunction<-function(distance){
  if(distance > 50e3)
    return(0)
  return(1/50e3*(50e3-distance))
}

scoreFunction<-function(R2, startPos, endPosVector){
  vec <- rep(0, length(R2))
  M <- 0.6
  vec[1] = scaleFunction(endPosVector[1] - startPos) * R2[1] - M
  if(length(R2) == 1)
    return(vec)
  
  for(i in 2:length(R2)){
    vec[i] = vec[i-1] + (scaleFunction(endPosVector[i] - endPosVector[i-1]) * R2[i] - M)
  }
  return(vec)
}

kernel<-function(vector, bandwidth = 5){
  normaliser<-2*sum(3/4*(1-((0:(bandwidth-1))/bandwidth)^2))
  
  ret<-rep(0,length(vector))
  for(i in (bandwidth+2):(length(vector)-bandwidth-1)){
    
    curSum = 0
    offset = bandwidth
    for(j in (i-(bandwidth+1)):(i-1)){
      curSum = curSum + vector[j]*(3/4*(1-((offset-1)/bandwidth)^2))/normaliser
      offset = offset - 1;
    }
    curSum = curSum + vector[i]
    offset = 1;
    for(j in (i+1):(i+(bandwidth+1))){
      curSum = curSum + vector[j]*(3/4*(1-((offset-1)/bandwidth)^2))/normaliser
      offset = offset + 1
    }
    #cat(i,"/",length(ret),": ", ret[i])
    #cat(curSum)
    
    ret[i] = curSum
  }
  return(ret)
}

plotRecomb<-function(data,chr,from,to,...){
  dat<-data[data$V2==chr,]
  plot(-1,-1,xlim=c(from,to),ylim=c(0,100), xaxs="i",yaxs="i", ...)
  for(i in 1:nrow(dat)){
    rect(dat$V3[i],0,dat$V4[i],dat$V5[i],col = c("black","grey")[as.factor(dat$gender)])
  }
}