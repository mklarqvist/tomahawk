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
mat2<-mat/1000 # Truncate at 1000
mat2[mat2>1]<-1 # Everything over 1 squash to 1
image(as.matrix(mat2),useRaster = T,axes=F)
axis(1, at = seq(0, 62943450, by = 10e6)/62943450, labels = seq(0, 62943450, by = 10e6)/1e6, las=2)
axis(2, at = seq(0, 62943450, by = 10e6)/62943450, labels = seq(0, 62943450, by = 10e6)/1e6, las=2)
