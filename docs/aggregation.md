# Aggregation in Tomahawk
## Motivation
[Tomahawk](https://github.com/mklarqvist/tomahawk) can output many millions to many hundreds of millions to billions of output linkage disequilibrium (LD) associations generated from many millions of input SNVs. Take for example a small chromosome like `chr20` using data from the 1000 Genomes Project Phase 3 (1KGP3). This data comprises of 2,504 samples with 1,733,484 diploid SNVs. Assuming we can plot the LD data for a pair of SNVs in a single pixel, a monitor would have to have the dimensions 400 x 400 meters to display this data*! Not only would the monitor be huge, the memory requirement for plotting this image would be around 400 GB! Here we describe methods to overcome these obstacles.

## Aggregation
Aggregation by summation

| | | | |
|----|----|----|----|
| 1  | 2  | 3  | 4  |
| 5  | 6  | 7  | 8  |
| 9  | 10 | 11 | 12 |
| 13 | 14 | 15 | 16 |

Can be reduced to

|    |    |
|----|----|
| 14 | 22 |
| 46 | 54 |


\* Assuming a 1920 x 1080 pixel resolution and 20" monitor as reference