# Job balancing in Tomahawk

## Motivation
[Tomahawk](https://github.com/mklarqvist/tomahawk) stores variants in non-overlapping blocks. These blocks has to be compared either to themselves (diagonal) or against each other (square) in order to compute linkage disequilibrium (LD) for a set of variants. First we describe the rationale for pre-loading data into memory and then how to partition this subset in the most efficient way.

## Pre-loading data
Without losing generality, consider the situation where you have a set of three blocks {1, 2, 3} that you want to compare pairwise. We describe the steps of a simple iterative algorithm to compare them below:
* Load blocks 1 and 2 into memory and compare {1, 2}
* Release block 2 from memory
* Load block 3 and compare {1, 3}
* Release block 1 from memory
* Load block 2 from memory and compare {2, 3}

Notice that we have now released and loaded the data for block 2 twice. This undesired memory- and IO-overhead grows square to the number of blocks, `O( (N-1) * (N-1) )` where `N` is the number of blocks. 
* If you have 10 blocks, there will be 9 overhead loads for each block for a total of 81 excess loads.
* If you have 100 blocks, there will be 99 overhead loads for each block for a total of 9,801 excess loads
* If you have 1000 blocks, there will be 999 overhead loads for each block for a total of 998,001 excess loads

For example, using standard import parameters, chromosome 20 for the 1000 Genomes Project data has 1,696 blocks. Without addressing this problem, we would have an excess of 2,873,025 overhead loads. Because of this exorbant cost it is very desirable to load all the data of interest into memory just once. However, if done without careful effort to partition data of interest into subproblems then memory will become limiting very quickly. We describe how we address this problem in Tomahawk below.

As an additional footnote, it is worthwile to mention that pre-loading data will spare the file-system on computer farms. This is generally always the most rate-limiting resource available.

## Memory-sparing job-loading
[Tomahawk](https://github.com/mklarqvist/tomahawk) splits large problems into multiple psuedo-balanced sub-problems in a memory-aware fashion using a tiling approach. 
Without losing generality, consider the situation where you want to calculate linkage disequilibrium (LD) for all variants pairwise. As a consequence, any given locus `v`
will be compared to every other loci, `V`. The simplest, naïve, way to parallelize this involves giving subproblem `j` out of `J` a list 
of loci from `[i*j, (i+1)*j]` for all subproblems, where `i = V/J`. This approach would invariantly require all data in memory irrespective of the slice-size. We can examplify this by drawing a square of four loci and divide the problem into three: {(1), (2), (3,4)}. Highlighted in bold are the sites addressed in subproblem 1 out of 3. Note that the lower triangular is not actually computed in practice and only shown here for visual clarity.

| Loci   | 1   | 2   | 3   | 4   |
|---|-----|-----|-----|-----|
| 1 | **1,1** | **2,1** | **3,1** | **4,1** |
| 2 | 1,2 | 2,2 | 3,2 | 4,2 |
| 3 | 1,3 | 2,3 | 3,3 | 4,3 |
| 4 | 1,4 | 2,4 | 3,4 | 4,4 |

In this approach, we need to have the data for variants {1, 2, 3, 4} when computing pairwise LD for locus {1}. In addition to the raw RLE objects, we would potentially have to allocate additional memory for the bit-vectors and masks. For most large datasets, this encompasses many tens of gigabytes. Not only does this require a large amount of memory but will also result in slower compute time because of the poor spatial locality of the data in memory.  

In order to overcome this restrictive boundary, we describe a slightly more complex load-balancing solution that revolves around the partitioning of the `V^2` problem space into a square grid of subproblems. This will always be possible in the complete case as both dimensions are equal (`V` and `V`). In this case, each subproblem gets a list of loci from `[i*j, (i+1)*j]` in the imaginary x-dimension and `[k*j, (k+1)*j]` in the y-dimension. Slicing the data in two dimensions enables us to reduce memory usage down from `O(V)` to `O(|k| + |l|)` in the worst case. We can examplify this approach using the same example from above

|   | 1   | 2   | 3   | 4   |
|---|-----|-----|-----|-----|
| 1 | **1,1** | **2,1** | 3,1 | 4,1 |
| 2 | **1,2** | **2,2** | 3,2 | 4,2 |
| 3 | 1,3 | 2,3 | 3,3 | 4,3 |
| 4 | 1,4 | 2,4 | 3,4 | 4,4 |

In this approach, we only need to have data for {1,2} when computing LD for {1,2}. This memory-sparing approach dramatically reduces the memory-requirements per job on most datasets. For this approach to generally work out, we need to choose a slice-size such that its upper triangular plus diagonal equals the number of subproblems we want to solve. Here's the sequence of the first 50 valid slice-sizes:  
1, 3, 6, 10, 15, 21, 28, 36, 45, 55, 66, 78, 91, 105, 120, 136, 153, 171, 190, 210, 231, 253, 276, 300, 325, 351, 378, 406, 435, 465, 496, 528, 561, 595, 630, 666, 703, 741, 780, 820, 861, 903, 946, 990, 1035, 1081, 1128, 1176, 1225, 1275  

You can generate this list in `R`:
```R
choose(1:50, 2) + 1:50 # triangular + diagonal
```
and make note that this is of course equivalent to
```R
( (1:50)^2 + (1:50) ) / 2 # (square + diagonal) / 2
```

The total number of subproblems to solve (`-c `) has to be a member of this function.

## Practical difference
We can demonstrate the efficiency of our grid-partitioning method even on small chromosomes like `chr20` using data from the 1000 Genomes Project Phase 3 (1KGP3). This data comprises of 2,504 samples with 1,733,484 diploid SNVs. As expected, when partitioning the problem as described we see a propotional decrease in memory allocation with the number of subproblems.

| Approach   | Memory  |
|------------|---------|
| Naïve      | 5 GB    |
| Tiling-10  | 1312 MB |
| Tiling-45  | 608 MB  |
| Tiling-105 | 398 MB  |

