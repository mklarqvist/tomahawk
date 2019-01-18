# Introduction

The primary goal of `tomahawk` is to efficiently and conveniently compute
linkage-disequilibrium for a single site vs its local neighbourhood, a sliding
window across an interval, or all-vs-all pairwise genome-wide. Tomahawk is a C++
library with a standard CLI divided into various subroutines, as is the standard
for bioinformatics tools.

*  `tomahawk` is *much* more efficient than existing solutions, both in terms of
   memory usage and simulation time. `tomahawk` can easily compute genome-wide
   LD for cohorts of millions of samples.

*  `tomahawk` is designed as a [C++ API](google.com) to simplify LD-based
    workflow: either directly by using the C++ API or using any of the available
    bindings. Currently there are [R bindings](www.google.com) and [Python3
    bindings](www.google.com).

*  `tomahawk` has its own human-readble and binary interchange format for
    expedient anaylsis of the output data. Other text-based systems are
    *extremely* inefficient and do not support basic operations such as
    subsetting, searching, summarizing, and visualizing.

*   `tomahawk` accepts any valid htslib-compatible input variant call format
    file for import into the internal tomahawk (`twk`) file format.