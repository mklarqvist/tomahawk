# CLI: `aggregate`
Short

## Synopsis
Aggregate TWO data into a rasterized matrix of size [x,y] for plotting.

## Examples
```bash
tomahawk aggregate -f r2 -r mean -c 50 -i input.twk
```

## Options
```bash
Options:
  -i   FILE   input TWO file (required)
  -x,y INT    number of X/Y-axis bins (default: 1000)
  -f   STRING aggregation function: can be one of (r2,r,d,dprime,dp,p,hets,alts,het,alt)(required)
  -r   STRING reduction function: can be one of (mean,count,n,min,max,sd)(required)
  -I   STRING filter interval <contig>:pos-pos (TWK/TWO) or linked interval <contig>:pos-pos,<contig>:pos-pos
  -c   INT    min cut-off value used in reduction function: value < c will be set to 0 (default: 5)
  -t   INT    number of parallel threads: each thread will use 40(x*y) bytes
```