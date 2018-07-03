# Aggregation in Tomahawk
## Motivation
[Tomahawk](https://github.com/mklarqvist/tomahawk) generally output many millions to many hundreds of millions to billions of output linkage disequilibrium (LD) associations generated from many millions of input SNVs. It is technically very challenging to visualize such large datasets&mdash;not only in terms of  

In order to get a scope at the scale of this problem, take for example a small chromosome like `chr20` using data from the 1000 Genomes Project Phase 3 (1KGP3). This data comprises of 2,504 samples with 1,733,484 diploid SNVs. Assuming we can plot the LD data for a pair of SNVs in a single pixel, a monitor would have to have the dimensions 400 x 400 meters to display this data*! Not only would the monitor have to be huge, the memory requirement for plotting this image would be around 400 GB! Here we describe methods to overcome these obstacles.

\* Assuming a 1920 x 1080 pixel resolution and 20" monitor as reference

## Existing solutions
There are several excellent solutions for rasterizing large datasets, including [Datashader](http://datashader.org/index.html). But packages like this requires us to leave the highly compressed binary representation in Tomahawk in order to transform `two` entries into a form understandable by these solutions. We have tried most of the popular applications that aggregate datasets and have found none that works in tiny memory and efficiently on our specific data.

## Aggregation
Aggregation, or rasterization, is the process of reducing larger datasets to smaller ones for the purposes of displaying more data than can fit on the screen at once. Tomahawk performs aggregation into regular grids by invoking some summary statistics function on your data in the given bins. At the moment, Tomahawk supports aggregation by
* summation
* mean
* minimum
* maximum
* standard deviation
* summation squared
* count

Without losing generality, imagine we start out with this 4x4 matrix of observations and we want to plot 4 pixels (2 x 2).

|    | C1 | C2 | C3 | C4 |
|----|----|----|----|----|
| **R1** | 1  | 2  | 3  | 4  |
| **R2** | 5  | 6  | 7  | 8  |
| **R3** | 9  | 10 | 11 | 12 |
| **R4** | 13 | 14 | 15 | 16 |

Aggregation by summation

|      | C1-2 | C3-4 |
|------|------|------|
| **R1-2** | 14   | 22   |
| **R3-4** | 46   | 54   |

Aggregation by mean

|      | C1-2 | C3-4 |
|------|------|------|
| **R1-2** | 3.5  | 5.5  |
| **R3-4** | 11.5 | 13.5 |

Aggregation by min

|      | C1-2 | C3-4 |
|------|------|------|
| **R1-2** | 1    | 3    |
| **R3-4** | 9    | 11   |

Aggregation by max

|      | C1-2 | C3-4 |
|------|------|------|
| **R1-2** | 6    | 8    |
| **R3-4** | 14   | 16   |

Aggregation by count

|      | C1-2 | C3-4 |
|------|------|------|
| **R1-2** | 4    | 4    |
| **R3-4** | 4    | 4    |

