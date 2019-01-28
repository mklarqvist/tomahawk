# Introduction

The primary goal of `tomahawk` is to efficiently and conveniently compute
linkage-disequilibrium for a single site vs its local neighbourhood, a sliding
window across an interval, or all-vs-all pairwise genome-wide. Tomahawk is a C++
library with a standard CLI divided into various subroutines, as is the standard
for most bioinformatics tools.

*  `tomahawk` is much more efficient than existing solutions, both in terms of
   memory usage and compute time. `tomahawk` can easily compute genome-wide
   LD for cohorts of millions of samples over chromosome-scaled regions.

*  `tomahawk` is designed as a [C++ API](api-documentation.md) to simplify LD-based
    workflow: either directly by using the C++ API or using any of the available
    language bindings. Currently there are [R bindings](r-tutorial.md) and [Python3
    bindings](python-tutorial.md).

*  `tomahawk` has its own human-readable (`.ld`) interchange format and a highly
    compressed binary (`.two`) format for expedient analysis of the generated
    output data. Other text-based systems are *extremely* inefficient and do not
    support basic operations such as subsetting, searching, summarizing, and
    visualizing.

*   `tomahawk` accepts any valid
    [htslib](https://github.com/samtools/htslib)-compatible input variant call
    format file for import into the internal binary `tomahawk` (`.twk`) file
    format.