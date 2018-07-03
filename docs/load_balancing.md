## Memory-sparing job-loading
[Tomahawk](https://github.com/mklarqvist/tomahawk) splits large problems into multiple psuedo-balanced sub-problems in a memory-aware fashion using a tiling approach. 
Without losing generality, consider the situation where you want to calculate LD for all variants pairwise. As a consequence, any given locus `v`
will be compared to every other loci, `V`. The simplest, naïve, way to parallelize this involves giving subproblem `j` out of `J` a list 
of loci from `[i*j, (i+1)*j]` for all subproblems. This approach would invariantly require all data in memory irrespective of the slice-size. We can examplify this by drawing a square of four loci and divide the problem into two. Highlighted in bold are the sites addressed in subproblem 1 out of 4.

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

In this approach, we only need to have data for {1,2} when computing LD for {1,2}. For this approach to work out, we need to choose a slice-size such that its upper triangular plus diagonal equals the number of subproblems we want to solve. Here's the sequence of the first 50 valid slice-sizes:  
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