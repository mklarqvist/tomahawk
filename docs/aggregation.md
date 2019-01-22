# Aggregating datasets

## Motivation

[Tomahawk](https://github.com/mklarqvist/tomahawk) generally output many
millions to many hundreds of millions, or even billions, of output linkage
disequilibrium (LD) associations generated from many millions of input SNVs. It
is technically very challenging to visualize such large datasets. Not only
because of hardware limitations such as loading all the data into memory, or
directly rendering billions of data points, but also because of more practical
considerations such as cramming such a vast number of data points into a finite
number of pixels would result in an absolute horrendous and uninformative image.

In order to get a scope of the scale this problem presents, take for example a
small chromosome like `chr20` with data from the 1000 Genomes Project Phase 3
(1KGP3). This data comprises of 1,733,484 diploid SNVs. Assuming we can plot the
LD data for a pair of SNVs in a single pixel, a monitor would have to have the
dimensions 400 x 400 meters to display this data*! Not only would the monitor
have to be huge, the memory requirement for plotting this image would be around
400 GB! Here we describe methods to overcome these obstacles.

\* Assuming a 1920 x 1080 pixel resolution and 20" monitor as reference

## Existing solutions
There are several existing solutions for aggregating large datasets, such as
[Datashader](http://datashader.org/) for Python users. But packages like this
requires us to leave the highly compressed internal binary representation of
Tomahawk in order to transform `two` records into a form understandable by these
frameworks. We have tried several of the most popular framework for aggregating
datasets and found none that works in reasonable memory and is sufficiently
efficient when applied to our specific use-case.

## Aggregation
Aggregation is the process of reducing larger datasets to smaller ones for the
purposes of displaying more data than can fit on the screen at once while
maintaining the primary features of the original dataset. Tomahawk performs
aggregation into regular grids (two-dimensional partitions) by applying summary
statistics function on data collected in the given bins. At the moment, Tomahawk
supports aggregation by

| Function           |  Action                                                              |
|--------------------|----------------------------------------------------------------------|
| Summation          | Sum total of the desired property                                    |
| Summation squared  | Sum total of squares of the desired property                         |
| Mean               | Mean of the desired property                                         |
| Standard deviation | Standard deviation of the desired property                           |
| Minimum            | Smallest value observed of the desired property                                                       |
| Maximum            | Largest value oserved of the desired property                                                        |
| Count              | Number of times a non-zero value is observed of the desired property |

Without losing generality, imagine we start out with this 4x4 matrix of
observations and we want to plot 4 pixels (2 x 2).

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

