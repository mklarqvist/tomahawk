color_range<-11
colors<-colorRampPalette(c("blue","red"))(color_range)
#colors[2:5]<-paste0(colors[2:5],c(20,40,60,80))
colors[1:2]<-"#FFFFFF"
mat<-read.delim("~/Downloads/1kgp3/1kgp3_chr20_matrix.txt",h=F)
# Linear transformation of data into percentile space for mapping the given
# colorkey to this set
#
# Transform matrix values into 1-percentile bins
dist<-table(cut(mat[mat>0],breaks=seq(0,max(mat[mat>0]),length.out = 101),include.lowest = T))
plot(cumsum(dist/sum(dist)),type="o",pch=20)
col_breaks<-rep(0,color_range-1)
for(i in 1:(color_range-1)){
  abline(v=which.max(cumsum(dist/sum(dist))>1 - 1/i))
  # Compute 10-percentile bins using the nearest-rank method
  col_breaks[i] = which.max(cumsum(dist/sum(dist))>1 - 1/i) / 100
}

image(as.matrix(mat),breaks = c(0, col_breaks, 1), col = colors,useRaster = T)

# For associative count matrix
mat<-read.delim("~/Downloads/1kgp3/1kgp3_chr20_matrix.txt",h=F)

jet.colors <-
  colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                     "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

# 
setwd("~/Desktop/1kgp3_aggregates/")
for(i in 1:22){
  tryCatch({
  mat<-read.delim(paste0("~/Downloads/1kgp3/chr",i,"_aggregate.out"),h=F,nrows = 4000)
  #mat2<-mat/round(mean(mat[mat>5])*5,-2) # Truncate at a count of 1000
  mat2<-mat/2000
  mat2[mat2>1]<-1 # Everything over 1 squash to 1
  
  
  dist<-table(cut(mat2[mat2>0],breaks=seq(0,max(mat2[mat2>0]),length.out = 101),include.lowest = T))
  #plot(cumsum(dist/sum(dist)),type="o",pch=20)
  col_breaks<-rep(0,color_range-1)
  for(j in 1:(color_range-1)){
    #abline(v=which.max(cumsum(dist/sum(dist)) > 1 - 1/j))
    # Compute 10-percentile bins using the nearest-rank method
    col_breaks[j] = which.max(cumsum(dist/sum(dist)) > 1 - 1/j) / 100
  }
  
  
  jpeg(paste0("1kgp3_chr",i,"_4k_aggregate_col4_linear_quantile_transform.jpeg"),width = 4000, height=4000)
  par(mar=c(0,0,0,0)) # set all margins to 0
  image(as.matrix(mat2),useRaster = T,axes=F, xaxt='n', yaxt='n',ann=FALSE, bty="n",col=viridis(11),breaks=c(0,col_breaks,1))
  dev.off()
  
  jpeg(paste0("1kgp3_chr",i,"_4k_aggregate_col4.jpeg"),width = 4000, height=4000)
  par(mar=c(0,0,0,0)) # set all margins to 0
  image(as.matrix(mat2),useRaster = T,axes=F, xaxt='n', yaxt='n',ann=FALSE, bty="n",col=viridis(11))
  dev.off()
  
  mat2<-mat/round(mean(mat[mat>5])*5,-2) # Truncate at a count of 1000
  mat2[mat2>1]<-1 # Everything over 1 squash to 1
  jpeg(paste0("1kgp3_chr",i,"_4k_aggregate_col4_scaled_local.jpeg"),width = 4000, height=4000)
  par(mar=c(0,0,0,0)) # set all margins to 0
  image(as.matrix(mat2),useRaster = T,axes=F, xaxt='n', yaxt='n',ann=FALSE, bty="n",col=viridis(11))
  dev.off()
  
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n"); })
}

setwd("~/Desktop/1kgp3_aggregates/subpopulations/")
files<-list.files("~/Downloads/1kgp3/subpopulations/",full.names = T)
pops<-strsplit("ACB,ASW,BEB,CDX,CEU,CHB,CHS,CLM,ESN,FIN,GBR,GIH,GWD,IBS,ITU,JPT,KHV,LWK,MSL,MXL,PEL,PJL,PUR,STU,TSI,YRI",",")[[1]]
for(i in 1:length(pops)){
  tryCatch({
    mat<-read.delim(files[i],h=F,nrows = 4000)
    #mat2<-mat/round(mean(mat[mat>5])*5,-2) # Truncate at a count of 1000
    mat2<-mat/2000
    mat2[mat2>1]<-1 # Everything over 1 squash to 1
    
    
    dist<-table(cut(mat2[mat2>0],breaks=seq(0,max(mat2[mat2>0]),length.out = 101),include.lowest = T))
    #plot(cumsum(dist/sum(dist)),type="o",pch=20)
    col_breaks<-rep(0,color_range-1)
    for(j in 1:(color_range-1)){
      #abline(v=which.max(cumsum(dist/sum(dist)) > 1 - 1/j))
      # Compute 10-percentile bins using the nearest-rank method
      col_breaks[j] = which.max(cumsum(dist/sum(dist)) > 1 - 1/j) / 100
    }
    
    
    jpeg(paste0("1kgp3_chr11_",pops[i],"_4k_aggregate_col4_linear_quantile_transform.jpeg"),width = 4000, height=4000)
    par(mar=c(0,0,0,0)) # set all margins to 0
    image(as.matrix(mat2),useRaster = T,axes=F, xaxt='n', yaxt='n',ann=FALSE, bty="n",col=viridis(11),breaks=c(0,col_breaks,1))
    dev.off()
    
    jpeg(paste0("1kgp3_chr11_",pops[i],"_4k_aggregate_col4.jpeg"),width = 4000, height=4000)
    par(mar=c(0,0,0,0)) # set all margins to 0
    image(as.matrix(mat2),useRaster = T,axes=F, xaxt='n', yaxt='n',ann=FALSE, bty="n",col=viridis(11))
    dev.off()
    
    mat2<-mat/round(mean(mat[mat>5])*5,-2) # Truncate at a count of 1000
    mat2[mat2>1]<-1 # Everything over 1 squash to 1
    jpeg(paste0("1kgp3_chr11_",pops[i],"_4k_aggregate_col4_scaled_local.jpeg"),width = 4000, height=4000)
    par(mar=c(0,0,0,0)) # set all margins to 0
    image(as.matrix(mat2),useRaster = T,axes=F, xaxt='n', yaxt='n',ann=FALSE, bty="n",col=viridis(11))
    dev.off()
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n"); })
}

axis(1, at = seq(0, 62943450, by = 10e6)/62943450, labels = seq(0, 62943450, by = 10e6)/1e6, las=2)
axis(2, at = seq(0, 62943450, by = 10e6)/62943450, labels = seq(0, 62943450, by = 10e6)/1e6, las=2)
