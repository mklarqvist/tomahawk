# R tutorial

## Plotting data with `rtomahawk`

Color schemes available
<img src="../images/rtwk_colors.jpeg">

### Plot square LD
```R
library(rtomahawk)
library(viridis)

# Load some local data into memory.
twk <- openTomahawkOutput("1kgp3_chr6.two")
y <- readRecords(twk,"6:5e6-10e6", really=TRUE)
# Plot two panels in the same figure.
# Top panel: opacity gradient from [0.1, 1.0] mapping to 
# R2 range [0.1, 0.1, 0.2, ..., 1.0]
# <COLOR> placeholder gets replaced with one of the color 
# schemes described below.
par(mfrow=c(2,1))
plotLD(y, ylim=c(5e6,8e6), xlim=c(5e6,8e6), colors=<COLOR>, bg=<COLOR>(11)[1])
plotLD(y, ylim=c(5e6,8e6), xlim=c(5e6,8e6), colors=<COLOR>, bg=<COLOR>(11)[1], opacity=FALSE)
```

| Color scheme |  Image |
|--------------|--------|
| `default` |  <img src="../images/twk_plotLD_default_quad.jpeg">      |
| `cividis`             | <img src="../images/twk_plotLD_cividis_quad.jpeg">       |
| `inferno`             | <img src="../images/twk_plotLD_inferno_quad.jpeg">       |
| `magma`             |<img src="../images/twk_plotLD_magma_quad.jpeg">        |
| `plasma`             |<img src="../images/twk_plotLD_plasma_quad.jpeg">        |
| `viridis`                |<img src="../images/twk_plotLD_viridis_quad.jpeg">        |

### Plot triangular LD
Default for range
```R
plotLDTriangular(y,colors=viridis(11),bg=viridis(11)[1])
```
<img src="../images/twk_plotLD_triangular.jpeg">

Truncate y-axis to show local neighbourhood
```R
plotLDTriangular(y,colors=viridis(11),bg=viridis(11)[1],ylim=c(0,300e3))
```

<img src="../images/twk_plotLD_triangular_truncate.jpeg">

Orientations
```R
par(mfrow=c(2,2))
plotLDTriangular(y,ylim=c(0,300e3),colors=viridis(11),bg=viridis(11)[1], orientation = 1)
plotLDTriangular(y,ylim=c(0,300e3),colors=viridis(11),bg=viridis(11)[1], orientation = 2)
plotLDTriangular(y,ylim=c(0,300e3),colors=viridis(11),bg=viridis(11)[1], orientation = 3)
plotLDTriangular(y,ylim=c(0,300e3),colors=viridis(11),bg=viridis(11)[1], orientation = 4)
```

<img src="../images/twk_plotLD_triangular_orientations.jpeg">
Without annotations
```R
plotLDTriangular(y,ylim=c(0,300e3),colors=viridis(11),bg=viridis(11)[1], orientation = 1, annotate = FALSE)
```
<img src="../images/twk_plotLD_triangular_orientations_no_annotation.jpeg">

Note that these plotting functions still respect the global `mar` (margin)
values and have white space around it. We can change the global `par` argument
to, for example, zero to remove these margins when annotation is disabled. 

<img src="../images/twk_plotLD_triangular_orientations_no_annotation_nomar.jpeg">

### LocusZoom
In this section we will produce LocusZoom-like plots using GWAS data from the UK
BioBank comprised exclusively of white British individuals. The Roslin Institute
at the University of Edinbrugh host a data browser of associatins called [Gene
Atlas](http://geneatlas.roslin.ed.ac.uk/). In the following examples we will
investigate the association of genotypes at chromosome 6 and diabetes mellitus
in this cohort. Data in its entirety can be [explored
further](http://geneatlas.roslin.ed.ac.uk/downloads/?traits=493) using the Gene
Atlas. To reproduce the results below download the 
[imputed](ttp://static.geneatlas.roslin.ed.ac.uk/gwas/allWhites/imputed/data.copy/imputed.allWhites.selfReported_n_1245.chr6.csv.gz) 
data for chromosome 6 and the associated 
[positional information](http://static.geneatlas.roslin.ed.ac.uk/gwas/allWhites/snps/extended/snps.imputed.chr6.csv.gz).

```R
# Downloaded data from:
# http://static.geneatlas.roslin.ed.ac.uk/gwas/allWhites/imputed/data.copy/imputed.allWhites.selfReported_n_1245.chr6.csv.gz
# http://static.geneatlas.roslin.ed.ac.uk/gwas/allWhites/snps/extended/snps.imputed.chr6.csv.gz
library(data.table) # For speedier reading of data.
# We zcat (uncompress gzipped archive) directly into fread.
x <- fread("zcat imputed.allWhites.selfReported_n_1245.chr6.csv.gz", sep=" ")
snp <- fread("zcat snps.imputed.chr6.csv.gz", sep=" ")
# Keep matching data available in both files. We neeed to maintain parity.
snp <- snp[match(x$SNP, snp$SNP),]
# Transform linear-scale P-value into -log10(P).
snp$p <- -log10(x$`PV-selfReported_n_1245`)

# Setup path to local Tomahawk file to compute LD against.
twk2<-new("twk")
twk2@file.path <- "1kgp3_chr6.twk"

# Load recombination data supplied with `rtomahawk`
data(gmap)
# Region of 1 Mb in either direction of target.
single <- plotLZ(twk2, "6:20694884", snp, gmap, window=1e6, minR2=0)
# Region of 50 Kb in either direction of target.
single <- plotLZ(twk2, "6:20694884", snp, gmap, window=50e3, minR2=0)
```
<img src="../images/twk_locuszoom.jpeg">
<img src="../images/twk_locuszoom_50k.jpeg">

These functions are extremely fast as the R-bindings use `.Call` commands to
communicate with the compiled C++ shared object:
```R
> system.time(plotLZ(twk2, "6:20694884", snp, gmap, window=1e6, minR2=0))
   user  system elapsed 
  1.405   0.001   1.407

> system.time(plotLZ(twk2, "6:20694884", snp, gmap, window=50e3, minR2=0))
   user  system elapsed 
  0.115   0.004   0.119
```
Almost all of this time is spent rendering symbols in the vectorized plot. The
same command using the CLI takes roughly 1/3rd of the time:
```bash
$ time twk scalc -i 1kgp3_chr6.twk -I 6:20694884 -w 1000000  > /dev/null
real  0m0.452s
user  0m0.937s
sys   0m0.061s
```

### Combining plots
`rtomahawk` renders plots using base-R.

```R
par(mfrow=c(2,1),=c(0,5,3,5))
single <- plotLZ(twk2, "6:20694884", snp, gmap, window=1e6, minR2=0, xlab="", xaxt="n")
twk<-openTomahawkOutput("~/Downloads/test_region.two")
y<-readRecords(twk,really=TRUE)
par(mar=c(5,5,0,5))
plotLDTriangular(y, ylim=c(0,500e3), xlim=c(19694884,21694884),colors=viridis(11),bg=viridis(11)[1], orientation = 2, cex=.25, main="")
```
<img src="../images/twk_locuszoom_combine.jpeg">

Add a gene track using data from `biomaRt` and drawn using `Sushi`:
```R
library(biomaRt)
mart <- useMart(host='http://grch37.ensembl.org',biomart="ensembl", dataset="hsapiens_gene_ensembl")
results <- getBM(attributes = c("ensembl_exon_id",
                    "exon_chrom_start","exon_chrom_end"),
                 filters = c("chromosome_name", "start", "end"), 
                 values=list(6, 20694884-1e6, 20694884+1e6),
                 mart=mart)

results2 <- getBM(attributes = c("ensembl_exon_id","hgnc_symbol", "chromosome_name", 
                                 "start_position", "end_position","strand"),
                 filters = c("chromosome_name", "start", "end"), 
                 values=list(6, 20694884-1e6, 20694884+1e6),
                 mart=mart)

results<-cbind(results,results2[,-1,drop=F])
results<-results[results$hgnc_symbol!="",]
results<-results[,c(2,3,4,5,8)]
names(results) <- c("start","stop","gene","chrom","strand")
results$score = "."
results<-results[,c("chrom","start","stop","gene","score", "strand")]

library(Sushi)
chrom = 6
chromstart = 20694884-1e6
chromend = 20694884+1e6
pg = plotGenes(geneinfo=results,chrom=chrom,chromstart=chromstart,chromend=chromend)

layout(matrix(1:3,ncol=1),heights = c(2,2,1))
par(mar=c(0,5,3,5))
single <- plotLZ(twk2, "6:20694884", snp, gmap, window=1e6, minR2=0, xlab="", xaxt="n")
twk<-openTomahawkOutput("~/Downloads/test_region.two")
y<-readRecords(twk,really=TRUE)
par(mar=c(5,5,0,5))
plotLDTriangular(y, ylim=c(0,500e3), xlim=c(19694884,21694884),colors=viridis(11),bg=viridis(11)[1], orientation = 2, cex=.25, main="")
par(mar=c(2,5,0,5))
pg = plotGenes(geneinfo=results,chrom=chrom,chromstart=chromstart,chromend=chromend,labeloffset=.5,fontsize=1,arrowlength = 0.025)
abline(v=20694884,lty="dashed",col="grey")
rtomahawk:::addGenomicAxis(c(chromstart,chromend),at = 1, las = 1, F)
```
<img src="../images/twk_locuszoom_combine_genes.jpeg">
<img src="../images/twk_locuszoom_combine_genes_zoom.jpeg"  >

Zoom into local region (100kb flanking region):
<img src="../images/twk_locuszoom_combine_100k.jpeg">

`rtomahawk` and `tomahawk` computes millions of LD associations and plots millions of data points in seconds on a single thread
```R
> system.time(f())
   user  system elapsed 
  6.336   0.260   6.233 
```

Controlling graphics
<img src="../images/twk_locuszoom_pch.jpeg">

### Aggregation

Quantile-normalized (left) or linear range (right)
```R
twk<-openTomahawkOutput("example.two")
x<-aggregate(twk,aggregation="r2", reduction="count",xbins=1000, ybins=1000,minCount=50, verbose=T, threads=8)
par(mar=c(5,5,5,8), mfrow=c(1,2))
plotAggregation(x, normalize = TRUE)
plotAggregation(x, normalize = FALSE)
```
<img src="../images/twk_aggregate_r2.jpeg">