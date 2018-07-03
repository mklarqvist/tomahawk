# Aggregation in Tomahawk
## Motivation
[Tomahawk](https://github.com/mklarqvist/tomahawk) can output many millions to many hundreds of millions to billions of output linkage disequilibrium (LD) associations generated from many millions of input SNVs. Take for example a small chromosome like `chr20` using data from the 1000 Genomes Project Phase 3 (1KGP3). This data comprises of 2,504 samples with 1,733,484 diploid SNVs. Assuming we can plot the LD data for a pair of SNVs in a single pixel, a monitor would have to have the dimensions 400 x 400 meters to display this data*! Not only would the monitor have to be huge, the memory requirement for plotting this image would be around 400 GB! Here we describe methods to overcome these obstacles.

## Aggregation
Aggregation is the process of reducing larger datasets to smaller ones for the purposes of displaying more data than can fit on the screen at once. While there are a multitude of ways to aggregate large datasets, we make significant use of summation and prioritization for numerical and categorical data, respectively. The following sections will provide examples for how we use aggregation to reduce the size of larger datasets.

Imagine we start out with this 4x4 matrix of observations

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


\* Assuming a 1920 x 1080 pixel resolution and 20" monitor as reference