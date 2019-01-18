# Algorithmic overview
Vectors of genotypes are represented as fixed-width run-length encoded (RLE)
objects. This encoding scheme is generally superior to dynamic-width encoding
approaches in terms of iteration speed (as no data processing is required) but
inferior in terms of compressibility (as bits are wasted). The word-width of the
RLE entries is fixed across a file and is determined contextually given the
total number of samples. 

We describe three efficient algorithms to calculate genome-wide linkage
disequilibrium for all pairwise alleles/genotypes in large-scale cohorts. The
algorithms exploit different concepts: 1) low genetic diversity and 2) large
memory registers on modern processors. 

1. The first algorithm directly compares fixed-width compressed RLE entries from
   two vectors in worst-case O(|RLE_A| + |RLE_B| + 1)-time.
2. The second transforms compressed RLE entries to uncompressed k-bit-vectors
   and use machine-optimized SIMD-instructions to horizontally compare two such
   bit-vectors in worst-case O(N/W)-time. This algorithm also exploits the
   relatively low genetic diversity within species using implicit heuristics. 
3. The third algorithm computes summary statistics only by maintaining a
   positional index of non-reference alleles and associated uncompressed
   1-bit-vectors for each genotypic vector in guaranteed
   O(min(|NON_REF_A|,|NON_REF_B))-time. The 1-bit vectors in this algorithm is
   different compared to the ones used in algorithm 2.
