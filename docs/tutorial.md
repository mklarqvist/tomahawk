# Getting started with `tomahawk`

This is an introductory tutorial for using Tomahawk. It will cover:

* Importing into Tomahawk
* Computing linkage-disequilibrium
* Aggregating and visualizing datasets
* Datasets

## Usage instructions
The CLI of Tomahawk comprises of several distinct subroutines (listed below).
Executing `tomahawk` gives a list of commands with brief descriptions and
`tomahawk <command>` gives detailed details for that command.

All primary Tomahawk commands operate on the binary Tomahawk `twk` and Tomahawk
output `two` file format. Interconversions between `twk` and `vcf`/`bcf` is
supported through the commands `import` for `vcf`/`bcf`->`twk` and `view` for
`twk`->`vcf`. Linkage disequilibrium data is written out in binary `two` format
or human-readable `ld` format.

| Command        | Description                                                 |
|----------------|-------------------------------------------------------------|
| [`aggregate`](cli/cli-aggregate)| data rasterization framework for `TWO` files                |
| [`calc`](cli/cli-calc)          | calculate linkage disequilibrium                            |
| [`scalc`](cli/cli-scalc)        | calculate linkage disequilibrium for a single site          |
| [`concat`](cli/cli-concat)      | concatenate `TWO` files from the same set of samples        |
| [`import`](cli/cli-import)      | import `VCF`/`VCF.gz`/`BCF` to `TWK`                        |
| [`sort`](cli/cli-sort)          | sort `TWO` file                                             |
| [`view`](cli/cli-view)          | `TWO`-&gt;`LD`/`TWO` view, `TWO` subset and filter          |
| [`haplotype`](cli/cli-haplotype)| extract per-sample haplotype strings in `FASTA`/binary format |
| [`relationship`](cli/cli-relationship) | compute marker-based pair-wise sample relationship matrices |
| [`decay`](cli/cli-decay)        | compute LD-decay over distance                              |
| [`prune`](cli/cli-prune)        | perform graph-based LD-pruning of variant sites             |

## Importing into Tomahawk
By design Tomahawk only operates on diploid and bi-allelic SNVs and as such
filters out indels and complex variants. Tomahawk does not support mixed phasing
of genotypes in the same variant (e.g. `0|0`, `0/1`). If mixed phasing is found
for a record, all genotypes for that site are converted to unphased genotypes.
This is a conscious design choice as this will internally invoke the correct
algorithm to use for mixed-phase cases.  

Importing standard files to Tomahawk involes using the `import` command. The
following command imports a `bcf` file and outputs `file.twk` while filtering
out variants with >20% missingness and sites that deviate from [Hardy-Weinberg](https://en.wikipedia.org/wiki/Hardy%E2%80%93Weinberg_principle)
equilibrium with a probability < 0.001.

```bash
tomahawk import -i file.bcf -o file -m 0.2 -H 1e-3
```

```bash
$ md5sum 1kgp3_chr6.bcf
8c05554ebe1a51e99be2471c96fad4d9  1kgp3_chr6.bcf
```

```bash
$ bcftools view 1kgp3_chr6.bcf -HG | wc -l
5024119
```

```bash
$ time tomahawk import -i 1kgp3_chr6.bcf -o 1kgp3_chr6
```

!!! Note "Auto-completion of file extensions"
    
    You do not have to append the `.twk` suffix to the output name as Tomahawk
    will automatically add this if missing. In the example above `"1kgp3_chr6"`
    will be converted to `"1kgp3_chr6.twk"` automatically. This is true for 
    most Tomahawk commands when using the CLI.

```text hl_lines="20 22 23 29 30"
Program:   tomahawk-264d039a-dirty (Tools for computing, querying and storing LD data)
Libraries: tomahawk-0.7.0; ZSTD-1.3.8; htslib 1.9
Contact: Marcus D. R. Klarqvist <mk819@cam.ac.uk>
Documentation: https://github.com/mklarqvist/tomahawk
License: MIT
----------
[2019-01-21 11:58:15,692][LOG] Calling import...
[2019-01-21 11:58:15,692][LOG][READER] Opening 1kgp3_chr6.bcf...
[2019-01-21 11:58:15,695][LOG][VCF] Constructing lookup table for 86 contigs...
[2019-01-21 11:58:15,695][LOG][VCF] Samples: 2,504...
[2019-01-21 11:58:15,695][LOG][WRITER] Opening 1kgp3_chr6.twk...
[2019-01-21 11:58:39,307][LOG] Duplicate site dropped: 6:18233985
[2019-01-21 11:59:26,749][LOG] Duplicate site dropped: 6:55137646
[2019-01-21 11:59:39,315][LOG] Duplicate site dropped: 6:67839893
[2019-01-21 11:59:47,176][LOG] Duplicate site dropped: 6:74373442
[2019-01-21 11:59:51,440][LOG] Duplicate site dropped: 6:77843171
[2019-01-21 12:00:41,784][LOG] Duplicate site dropped: 6:121316830
[2019-01-21 12:01:12,830][LOG] Duplicate site dropped: 6:148573620
[2019-01-21 12:01:16,557][LOG] Duplicate site dropped: 6:151397786
[2019-01-21 12:01:42,644][LOG] Wrote: 4,784,608 variants to 9,570 blocks...
[2019-01-21 12:01:42,644][LOG] Finished: 03m26,952s
[2019-01-21 12:01:42,644][LOG] Filtered out 239,511 sites (4.76722%):
[2019-01-21 12:01:42,645][LOG]    Invariant: 15,485 (0.308213%)
[2019-01-21 12:01:42,645][LOG]    Missing threshold: 0 (0%)
[2019-01-21 12:01:42,645][LOG]    Insufficient samples: 0 (0%)
[2019-01-21 12:01:42,645][LOG]    Mixed ploidy: 0 (0%)
[2019-01-21 12:01:42,645][LOG]    No genotypes: 0 (0%)
[2019-01-21 12:01:42,645][LOG]    No FORMAT: 0 (0%)
[2019-01-21 12:01:42,645][LOG]    Not biallelic: 26,277 (0.523017%)
[2019-01-21 12:01:42,645][LOG]    Not SNP: 194,566 (3.87264%)

```

## Computing linkage-disequilibrium

### All pairwise comparisons

In this first example, we will compute all-vs-all LD associations using the data we imported
in the previous section. To limit compute, we restrict our attention to a 1/45 section of the
data by passing the `-c` and `-C` job parameters. We will be using 8 threads (`-t`), but you
may need to modify this to match the hardware available on your host machine. This job involves
comparing a pair of variants >141 billion times and as such takes around 30 min to finish on
most machines.
```bash
$ tomahawk calc -pi 1kgp3_chr6.twk -o 1kgp3_chr6_1_45 -C 1 -c 45 -t 8
```
!!! Note "Valid job balancing partitions"

    When computing genome-wide LD the
    balancing requires that number of sub-problems (`-c`) must be a member of the function: 
    $$
    \binom{c}{2} + c = \frac{c^2 + c}{2}, c > 0
    $$
    This function is equivalent to the upper-triangular of a square
    (`c`-by-`c`) matrix plus the diagonal. Read more about [load
    partitioning](job-balancing.md) in Tomahawk.

    Here is the first 100 valid partition sizes:  
    `1`, `3`, `6`, `10`, `15`, `21`, `28`, `36`, `45`, `55`, `56`, `68`, `81`, `95`, `110`, `126`, `143`, `161`, `180`, `200`, `211`, `233`, `256`, `280`, `305`, `331`, `358`, `386`, `415`, `445`, `466`, `498`, `531`, `565`, `600`, `636`, `673`, `711`, `750`, `790`, `821`, `863`, `906`, `950`, `995`, `1041`, `1088`, `1136`, `1185`, `1235`, `1276`, `1328`, `1381`, `1435`, `1490`, `1546`, `1603`, `1661`, `1720`, `1780`, `1831`, `1893`, `1956`, `2020`, `2085`, `2151`, `2218`, `2286`, `2355`, `2425`, `2486`, `2558`, `2631`, `2705`, `2780`, `2856`, `2933`, `3011`, `3090`, `3170`, `3241`, `3323`, `3406`, `3490`, `3575`, `3661`, `3748`, `3836`, `3925`, `4015`, `4096`, `4188`, `4281`, `4375`, `4470`, `4566`, `4663`, `4761`, `4860`, `4960`

!!! Warning "Massive workload"
    
    By default, Tomahawk will attempt to compute the genome-wide linkage-disequilibrium
    for all the provided variants and samples. This generally require huge amount
    of pairwise variant comparisons. Depending on your interests, you may not want to
    undertake this large compute and need to parameterize for your particular use-case.

!!! Hint "Different load balancing"
    
    As a consequence of the job balancing approach in Tomahawk, the diagonal jobs
    will always involve less variant comparisons (\\( \binom{n}{2} \\)) compared to off-diagonal
    jobs (\\(n \times m\\)), where \\(n\\) and \\(m\\) are the number of records in either block. 
    Keep this in mind if you are choosing a partition size given
    the run-times of a given job. For example, examine the total run-time of the second
    job (`-C 2`) &mdash; the first off-diagonal job &mdash; as a proxy for the expected average run-time per job.

!!! Note "Output ordering"
    
    By design, Tomahawk computes pairwise associations using an [out-of-order execution](https://simple.wikipedia.org/wiki/Out-of-order_execution) 
    paradigm and further permutes the input ordering by using [cache-blocking](https://en.wikipedia.org/wiki/Loop_nest_optimization).
    Additionally, each thread/slave transiently store a number of computed associations in a
    private buffer that will be flushed upon reaching some frequency threshold. Because of these
    technicalities, Tomahawk will neither produce ordered output nor will the permutation of the
    output order be identical between runs on identical data. This has no practical consequence
    to most downstream applications in Tomahawk with the exception of query speeds. To address
    this limitation we have introduced a powerful sorting paradigm that will be discussed below.

```text hl_lines="9 10 16 18"
Program:   tomahawk-264d039a-dirty (Tools for computing, querying and storing LD data)
Libraries: tomahawk-0.7.0; ZSTD-1.3.8; htslib 1.9
Contact: Marcus D. R. Klarqvist <mk819@cam.ac.uk>
Documentation: https://github.com/mklarqvist/tomahawk
License: MIT
----------
[2019-01-21 12:13:50,210][LOG] Calling calc...
[2019-01-21 12:13:50,211][LOG][READER] Opening 1kgp3_chr6.twk...
[2019-01-21 12:13:50,212][LOG] Samples: 2,504...
[2019-01-21 12:13:50,212][LOG][BALANCING] Using ranges [0-1063,0-1063] in square mode...
[2019-01-21 12:13:50,212][LOG] Allocating 1,063 blocks...
[2019-01-21 12:13:50,213][LOG] Running in standard mode. Pre-computing data...
[2019-01-21 12:13:50,213][LOG][SIMD] Vectorized instructions available: SSE4...
[2019-01-21 12:13:50,213][LOG] Constructing list, vector, RLE...
[2019-01-21 12:13:50,213][LOG][THREAD] Unpacking using 7 threads: ....... Done! 01,291s
[2019-01-21 12:13:51,504][LOG] 531,500 variants from 1,063 blocks...
[2019-01-21 12:13:51,504][LOG][PARAMS] square=TRUE,window=FALSE,low_memory=FALSE,bitmaps=FALSE,single=FALSE,force_phased=TRUE,force_unphased=FALSE,compression_level=1,block_size=500,output_block_size=10000,l_surrounding=500000,minP=1.000000,minR2=0.100000,maxR2=100.000000,minDprime=0.000000,maxDprime=100.000000,n_chunks=45,c_chunk=0,n_threads=8,ldd_type=3,cycle_threshold=0
[2019-01-21 12:13:51,504][LOG] Performing: 141,245,859,250 variant comparisons...
[2019-01-21 12:13:51,504][LOG][WRITER] Opening 1kgp3_chr6_1_45.two...
[2019-01-21 12:13:51,505][LOG][THREAD] Spawning 8 threads: ........
[2019-01-21 12:13:51,533][PROGRESS] Time elapsed       Variants           Genotypes         Output  Progress	Est. Time left
[2019-01-21 12:14:21,533][PROGRESS]      30,000s  3,685,746,500   9,229,109,236,000        988,858   2.60945%	18m39s
[2019-01-21 12:14:51,533][PROGRESS]   01m00,000s  6,425,868,750  16,090,375,350,000      1,689,422   4.54942%	20m58s
[truncated]
[2019-01-21 12:39:51,546][PROGRESS]   26m00,013s   140,131,382,750 350,888,982,406,000     48,635,114    99.211%	12s
[2019-01-21 12:40:04,317][PROGRESS] Finished in 26m12,784s. Variants: 141,245,859,250, genotypes: 353,679,631,562,000, output: 49,870,388
[2019-01-21 12:40:04,317][PROGRESS] 89,806,242 variants/s and 224,874,830,855 genotypes/s
[2019-01-21 12:40:04,321][LOG][PROGRESS] All done...26m12,815s!
```

This run generated 935817820 bytes (935.8 MB) of output data in binary format using the default
compression level (1):
```bash
$ ls -l 1kgp3_chr6_1_45.two
-rw-rw-r-- 1 mk819 mk819 935817820 Jan 21 12:40 1kgp3_chr6_1_45.two
```

If this data was stored in plain, human-readable, text format it would use 
over 45 GB:
```bash
$ tomahawk view -i 1kgp3_chr6_1_45.two | wc -c
4544702179
```

### Sliding window
If you are working with a species with well-know LD structure (such as humans) you can reduce the computational cost by limiting the search-space to a fixed-sized sliding window (`-w`). In window mode you are free to choose any arbitrary sub-problem (`-c`) size.

!!! Note "Window properties"
    
    By default, Tomahawk computes the linkage disequilibrium associations for all pairs
    of variants within `-w` bases of the target SNV on either direction. In order to
    maximize computational throughput, we made the design decision to *maintain* blocks
    of data that overlap to the target interval. This has the consequence that adjacent
    data that may **not** directly overlap with the target region will still be computed. The 
    technical reasoning for this is, in simplified terms, that the online repackaging of internal data
    blocks is more expensive than computing a small (relatively) number of non-desired (off-target) associations.

```bash
time tomahawk calc -pi 1kgp3_chr6.twk -o 1kgp3_chr6_4mb -w 4000000
```

```text
Program:   tomahawk-264d039a-dirty (Tools for computing, querying and storing LD data)
Libraries: tomahawk-0.7.0; ZSTD-1.3.8; htslib 1.9
Contact: Marcus D. R. Klarqvist <mk819@cam.ac.uk>
Documentation: https://github.com/mklarqvist/tomahawk
License: MIT
----------
[2019-01-21 14:10:33,616][LOG] Calling calc...
[2019-01-21 14:10:33,616][LOG][READER] Opening 1kgp3_chr6.twk...
[2019-01-21 14:10:33,619][LOG] Samples: 2,504...
[2019-01-21 14:10:33,619][LOG][BALANCING] Using ranges [0-9570,0-9570] in window mode...
[2019-01-21 14:10:33,619][LOG] Allocating 9,570 blocks...
[2019-01-21 14:10:33,620][LOG] Running in standard mode. Pre-computing data...
[2019-01-21 14:10:33,620][LOG][SIMD] Vectorized instructions available: SSE4...
[2019-01-21 14:10:33,620][LOG] Constructing list, vector, RLE...
[2019-01-21 14:10:33,620][LOG][THREAD] Unpacking using 7 threads: ....... Done! 15,774s
[2019-01-21 14:10:49,394][LOG] 4,784,608 variants from 9,570 blocks...
[2019-01-21 14:10:49,394][LOG][PARAMS] square=TRUE,window=TRUE,low_memory=FALSE,bitmaps=FALSE,single=FALSE,force_phased=TRUE,force_unphased=FALSE,compression_level=1,block_size=500,output_block_size=10000,window_size=4000000,l_surrounding=500000,minP=1.000000,minR2=0.100000,maxR2=100.000000,minDprime=0.000000,maxDprime=100.000000,n_chunks=1,c_chunk=0,n_threads=8,ldd_type=3,cycle_threshold=0
[2019-01-21 14:10:49,394][LOG] Performing: 11,446,234,464,528 variant comparisons...
[2019-01-21 14:10:49,394][LOG][WRITER] Opening 1kgp3_chr6_4mb.two...
[2019-01-21 14:10:49,402][LOG][THREAD] Spawning 8 threads: ........
[2019-01-21 14:10:49,424][PROGRESS] Time elapsed       Variants           Genotypes         Output  Progress	Est. Time left
[2019-01-21 14:11:19,424][PROGRESS]      30,000s  2,311,240,500   5,787,346,212,000        853,368         0	0
[2019-01-21 14:11:49,424][PROGRESS]   01m00,000s  4,562,731,500  11,425,079,676,000      1,717,956         0	0
[truncated]
[2019-01-21 16:12:49,492][PROGRESS] 02h02m00,068s   527,889,726,500  1,321,835,875,156,000    470,354,846         0	0
[2019-01-21 16:13:17,236][PROGRESS] Finished in 02h02m27,812s. Variants: 529,807,522,528, genotypes: 1,326,638,036,410,112, output: 473,514,018
[2019-01-21 16:13:17,237][PROGRESS] 72,104,114 variants/s and 180,548,701,880 genotypes/s
[2019-01-21 16:13:17,257][LOG][PROGRESS] All done...02h02m27,854s!
```

## Concatenating multiple archives

On of the immediate downsides of partitioning compute into multiple non-overlapping 
sub-problems is that we will generate a large number of independent files that must
be merged prior to downstream analysis. Fortunatly, this is a trivial operation in
Tomahawk and involves the `concat` command.

First lets compute some data using the first 3 out of 990 partitions of the dataset
used above.
```bash
for i in {1..3}; do time tomahawk calc -pi 1kgp3_chr6.twk -c 990 -C $i -o part$i\_3.two; done
```

Next, since we only have three files, we can concatenate (merge) these files together
into a single archieve by passing each file name to the command:
```bash
$ time tomahawk concat -i part1_3.two -i part2_3.two -i part3_3.two -o part3_concat
[2019-01-22 10:55:47,799][LOG] All files are compatible. Beginning merging...
[2019-01-22 10:55:47,800][LOG] Appending part1_3.two... 397.953344 Mb/98.987848 Mb
[2019-01-22 10:55:47,889][LOG] Appending part2_3.two... 359.538884 Mb/54.763322 Mb
[2019-01-22 10:55:47,938][LOG] Appending part3_3.two... 346.075500 Mb/52.072741 Mb
[2019-01-22 10:55:47,988][LOG] Finished. Added 3 files...
[2019-01-22 10:55:47,988][LOG] Total size: Uncompressed = 1.103568 Gb and compressed = 205.823911 Mb

real	0m0.226s
user	0m0.036s
sys	    0m0.190s
```

Passing every single input file in the command line does not become feasble when we have
many hundreds to thousands of files. To address this, it is possible to first store the
target file paths in a text file and then pass that to Tomahawk. Using the data from above,
we save these file paths to the file `part_file_list.txt`.
```bash
for i in {1..3}; do echo part$i\_3.two >> part_file_list.txt; done
```

Then we simply pass this list of file paths to Tomahawk using the `-I` argument:
```bash
$ time tomahawk concat -I part_file_list.txt -o part3_concat 
adding file=part1_3.two
adding file=part2_3.two
adding file=part3_3.two
[2019-01-22 10:57:43,227][LOG] All files are compatible. Beginning merging...
[2019-01-22 10:57:43,262][LOG] Appending part1_3.two... 397.953344 Mb/98.987848 Mb
[2019-01-22 10:57:43,350][LOG] Appending part2_3.two... 359.538884 Mb/54.763322 Mb
[2019-01-22 10:57:43,400][LOG] Appending part3_3.two... 346.075500 Mb/52.072741 Mb
[2019-01-22 10:57:43,447][LOG] Finished. Added 3 files...
[2019-01-22 10:57:43,447][LOG] Total size: Uncompressed = 1.103568 Gb and compressed = 205.823911 Mb

real	0m0.477s
user	0m0.040s
sys	    0m0.241s
```

!!! Danger "Possible duplications"
    
    The `concat` subroutine will dedupe input file paths but will **neiter check nor
    guarantee** that the actual input data is not duplicated. Therefore, it is possible to get
    duplicates in your data if you are not careful. These duplicate entries could corrupt
    any downstream insights!

## Sorting output files
All subroutines in Tomahawk will work on unsorted output files. However, sorted
files are *much* faster to query and will result in considerable time savings.
This is especially true for very large files. Sorting files in Tomahawk is
trivial and involves calling the `sort` command with the required fields input
(`-i`) and output (`-o`).

!!! Warning "Disk usage"
    
    Tomahawk can easily sort huge `.two` files with the help of [external merge
    sorting](https://en.wikipedia.org/wiki/External_sorting) algorithms followed by
    a single [k-way merge](https://en.wikipedia.org/wiki/K-way_merge_algorithm),
    memory assisted operation. The external merge step require additional disk space
    approximately equal to the size of the input dataset. Then, the merge operation
    will write a, generally, smaller output file. On average, we can approxiate that
    the sort operation require an additional two times the input file size on disk.
    After the sorting procedure is completed, the intermediate files are removed.

!!! Warning "Temporary files"
    
    By default, Tomahawk will generate temporary files in the same directory as the
    input file. If this file path is undesired then you have to modify the appropriate
    parameters.

In this example we will sort the almost 500 million (473,514,826) records
computed from the sliding window example above. This archive corresponds to >50
Gb of binary data.
```bash
$ time tomahawk sort -i 1kgp3_chr6_4mb.two -o 1kgp3_chr6_4mb_sorted
```

```text
Program:   tomahawk-264d039a-dirty (Tools for computing, querying and storing LD data)
Libraries: tomahawk-0.7.0; ZSTD-1.3.8; htslib 1.9
Contact: Marcus D. R. Klarqvist <mk819@cam.ac.uk>
Documentation: https://github.com/mklarqvist/tomahawk
License: MIT
----------
[2019-01-22 11:20:43,978][LOG] Calling sort...
[2019-01-22 11:20:43,993][LOG] Blocks: 47,358
[2019-01-22 11:20:43,993][LOG] Uncompressed size: 50.192950 Gb
[2019-01-22 11:20:43,993][LOG] Sorting 473,514,826 records...
[2019-01-22 11:20:43,993][LOG][THREAD] Data/thread: 6.274119 Gb
[2019-01-22 11:20:43,993][PROGRESS] Time elapsed       Variants  Progress	Est. Time left
[2019-01-22 11:20:43,994][LOG][THREAD] Slave-0: range=0->5919/47358 and name 1kgp3_chr6_4mb_sorted_fileSucRkC.two
[2019-01-22 11:20:43,994][LOG][THREAD] Slave-1: range=5919->11838/47358 and name 1kgp3_chr6_4mb_sorted_filewPI10T.two
[2019-01-22 11:20:43,994][LOG][THREAD] Slave-2: range=11838->17757/47358 and name 1kgp3_chr6_4mb_sorted_fileo1zcHb.two
[2019-01-22 11:20:43,994][LOG][THREAD] Slave-3: range=17757->23676/47358 and name 1kgp3_chr6_4mb_sorted_fileBvRnnt.two
[2019-01-22 11:20:43,994][LOG][THREAD] Slave-4: range=23676->29595/47358 and name 1kgp3_chr6_4mb_sorted_filewbBz3K.two
[2019-01-22 11:20:43,994][LOG][THREAD] Slave-5: range=29595->35514/47358 and name 1kgp3_chr6_4mb_sorted_fileCrNLJ2.two
[2019-01-22 11:20:43,994][LOG][THREAD] Slave-6: range=35514->41433/47358 and name 1kgp3_chr6_4mb_sorted_filerCfYpk.two
[2019-01-22 11:20:43,994][LOG][THREAD] Slave-7: range=41433->47358/47358 and name 1kgp3_chr6_4mb_sorted_filexx9a6B.two
[2019-01-22 11:21:13,994][PROGRESS]      30,000s     70,594,275   14.9086%	2m51s
[2019-01-22 11:21:43,994][PROGRESS]   01m00,000s    146,732,835    30.988%	2m13s
[2019-01-22 11:22:13,994][PROGRESS]   01m30,000s    245,535,675   51.8539%	1m23s
[2019-01-22 11:22:43,995][PROGRESS]   02m00,001s    323,064,235   68.2268%	55s
[2019-01-22 11:23:13,995][PROGRESS]   02m30,001s    421,592,790   89.0348%	18s
[2019-01-22 11:23:24,951][LOG][THREAD] Finished: 1kgp3_chr6_4mb_sorted_filewbBz3K.two with 23676-29595. Sorted n=59190000 variants with size=1.176156 Gb
[2019-01-22 11:23:25,114][LOG][THREAD] Finished: 1kgp3_chr6_4mb_sorted_fileBvRnnt.two with 17757-23676. Sorted n=59190000 variants with size=1.286938 Gb
[2019-01-22 11:23:25,385][LOG][THREAD] Finished: 1kgp3_chr6_4mb_sorted_filerCfYpk.two with 35514-41433. Sorted n=59190000 variants with size=1.219126 Gb
[2019-01-22 11:23:25,710][LOG][THREAD] Finished: 1kgp3_chr6_4mb_sorted_fileCrNLJ2.two with 29595-35514. Sorted n=59190000 variants with size=1.161081 Gb
[2019-01-22 11:23:26,384][LOG][THREAD] Finished: 1kgp3_chr6_4mb_sorted_filexx9a6B.two with 41433-47358. Sorted n=59184826 variants with size=1.215231 Gb
[2019-01-22 11:23:26,425][LOG][THREAD] Finished: 1kgp3_chr6_4mb_sorted_fileSucRkC.two with 0-5919. Sorted n=59190000 variants with size=1.242923 Gb
[2019-01-22 11:23:29,966][LOG][THREAD] Finished: 1kgp3_chr6_4mb_sorted_fileo1zcHb.two with 11838-17757. Sorted n=59190000 variants with size=1.686758 Gb
[2019-01-22 11:23:31,300][LOG][THREAD] Finished: 1kgp3_chr6_4mb_sorted_filewPI10T.two with 5919-11838. Sorted n=59190000 variants with size=2.013554 Gb
[2019-01-22 11:23:31,338][PROGRESS] 02m47,344s	473,514,826 (2,829,587 variants/s)
[2019-01-22 11:23:31,338][PROGRESS] Finished!
[2019-01-22 11:23:31,341][LOG] Spawning 112 queues with 2.380952 Mb each...
[2019-01-22 11:23:35,090][LOG][WRITER] Opening "1kgp3_chr6_4mb_sorted.two"...
[2019-01-22 11:23:35,091][PROGRESS] Time elapsed       Variants  Progress	Est. Time left
[2019-01-22 11:24:05,091][PROGRESS]      30,000s     39,757,023   8.39615%	5m27s
[2019-01-22 11:24:35,091][PROGRESS]   01m00,000s     76,990,000   16.2593%	5m9s
[2019-01-22 11:25:05,091][PROGRESS]   01m30,000s    107,840,000   22.7744%	5m5s
[2019-01-22 11:25:35,091][PROGRESS]   02m00,000s    133,420,000   28.1765%	5m5s
[2019-01-22 11:26:05,092][PROGRESS]   02m30,000s    170,215,235   35.9472%	4m27s
[2019-01-22 11:26:35,092][PROGRESS]   03m00,001s    208,070,000   43.9416%	3m49s
[2019-01-22 11:27:05,092][PROGRESS]   03m30,001s    247,010,000   52.1652%	3m12s
[2019-01-22 11:27:35,092][PROGRESS]   04m00,001s    287,690,000   60.7563%	2m35s
[2019-01-22 11:28:05,092][PROGRESS]   04m30,001s    325,480,000    68.737%	2m2s
[2019-01-22 11:28:35,093][PROGRESS]   05m00,001s    364,299,966   76.9353%	1m29s
[2019-01-22 11:29:05,093][PROGRESS]   05m30,002s    401,367,427   84.7634%	59s
[2019-01-22 11:29:35,093][PROGRESS]   06m00,002s    439,290,000   92.7722%	28s
[2019-01-22 11:30:02,180][PROGRESS] 06m27,089s	473,514,826 (1,223,269 variants/s)
[2019-01-22 11:30:02,180][PROGRESS] Finished!
[2019-01-22 11:30:02,194][LOG] Finished merging! Time: 06m27,103s
[2019-01-22 11:30:02,194][LOG] Deleting temp files...
[2019-01-22 11:30:02,194][LOG] Deleted 1kgp3_chr6_4mb_sorted_fileSucRkC.two
[2019-01-22 11:30:02,194][LOG] Deleted 1kgp3_chr6_4mb_sorted_filewPI10T.two
[2019-01-22 11:30:02,194][LOG] Deleted 1kgp3_chr6_4mb_sorted_fileo1zcHb.two
[2019-01-22 11:30:02,194][LOG] Deleted 1kgp3_chr6_4mb_sorted_fileBvRnnt.two
[2019-01-22 11:30:02,194][LOG] Deleted 1kgp3_chr6_4mb_sorted_filewbBz3K.two
[2019-01-22 11:30:02,194][LOG] Deleted 1kgp3_chr6_4mb_sorted_fileCrNLJ2.two
[2019-01-22 11:30:02,194][LOG] Deleted 1kgp3_chr6_4mb_sorted_filerCfYpk.two
[2019-01-22 11:30:02,194][LOG] Deleted 1kgp3_chr6_4mb_sorted_filexx9a6B.two
[2019-01-22 11:30:02,678][LOG] Finished!

real	9m18.793s
user	23m31.390s
sys	    0m53.197s

```

!!! success "Faster queries"
    
    Success: You can now use the newly created sorted archieve exactly as unsorted
    files but with faster query times.

## Format descriptions

!!! Note "Memory layout (technical)"

    A frequent design choice when designing a file format is how to align memory of records:
    either as lists of records (array-of-struct) versus column-orientated layouts (struct-of-array).
    The pivoted layout (struct-of-array) of `twk1_two_t` structs in Tomahawk will result *considerable*
    savings in disk space usage at the expense of querying speed of the resulting data. We decided to
    build the `two` format around the classical array-of-struct memory layout as we want to maximize
    the computability of the data. This is especially true in our application as we will never
    support individual columnar slicing and subset operations.

### `LD` format
Tomahawk can output binary `two` data in the human-readable `ld` format by invoking the `view` command. The primary output columns are described below:

| Column    | Description |
|----------|-------------|
| `FLAG`     | Bit-packed boolean flags (see below) |
| `CHROM_A`  | Chromosome for marker A            |
| `POS_A`    | Position for marker A            |
| `CHROM_B`  | Chromosome for marker B            |
| `POS_B`    | Position for marker B            |
| `REF_REF`  | Inner product of (0,0) haplotypes            |
| `REF_ALT`  | Inner product of (0,1) haplotypes            |
| `ALT_REF`  | Inner product of (1,0) haplotypes            |
| `ALT_ALT`  | Inner proruct of (1,1) haplotypes            |
| [`D`](https://en.wikipedia.org/wiki/Linkage_disequilibrium)        | Coefficient of linkage disequilibrium            |
| `DPrime`   | Normalized coefficient of linkage disequilibrium (scaled to [-1,1])            |
| [`R`](https://en.wikipedia.org/wiki/Pearson_correlation_coefficient)        | Pearson correlation coefficient            |
| `R2`       | Squared pearson correlation coefficient            |
| [`P`](https://en.wikipedia.org/wiki/Fisher's_exact_test)        | Fisher's exact test P-value of the 2x2 haplotype contigency table            |
| `ChiSqModel` | Chi-squared critical value of the 3x3 unphased table of the selected cubic root (&alpha;, &beta;, or &delta;)            |
| `ChiSqTable` | Chi-squared critical value of table (useful comparator when `P` = 0)            |

The 2x2 contingency table, or matrix, for the Fisher's exact test (`P`) for haplotypes look like this:

|                | REF-A | REF-B |
|----------------|---------------|-------------------|
| **REF-B** | A             | B                 |
| **ALT-B** | C             | D                 |

The 3x3 contigency table, or matrix, for the Chi-squared test for the unphased model looks like this:

|     | 0/0 | 0/1 | 1/1 |
|-----|-----|-----|-----|
| **0/0** | A   | B   | C   |
| **0/1** | D   | E   | F   |
| **1/1** | G   | H   | J   |

The `two` `FLAG` values are bit-packed booleans in a single integer field and describe a variety of states a pair of markers can be in.

| Bit position   | Numeric value | One-hot           | Description                                                                           |
|----------------|---------------|-------------------|---------------------------------------------------------------------------------------|
| 1              | 1             | 000000000000000**1** | Used phased math.                                                                     |
| 2              | 2             | 00000000000000**1**0 | Acceptor and donor variants are on the same contig.                                   |
| 3              | 4             | 0000000000000**1**00 | Acceptor and donor variants are far apart on the same contig.                         |
| 4              | 8             | 000000000000**1**000 | The output contingency matrix has at least one empty cell (referred to as complete).  |
| 5              | 16            | 00000000000**1**0000 | Output correlation coefficient is perfect (1.0).                                      |
| 6              | 32            | 0000000000**1**00000 | Output solution is one of >1 possible solutions. This only occurs for unphased pairs. |
| 7              | 64            | 000000000**1**000000 | Output data was generated in 'fast mode'.                                             |
| 8              | 128           | 00000000**1**0000000 | Output data is estimated from a subsampling of the total pool of genotypes.           |
| 9              | 256           | 0000000**1**00000000 | Donor vector has missing value(s).                                                    |
| 10             | 512           | 000000**1**000000000 | Acceptor vector has missing value(s).                                                 |
| 11             | 1024          | 00000**1**0000000000 | Donor vector has low allele count (<5).                                               |
| 12             | 2048          | 0000**1**00000000000 | Acceptor vector has low allele count (<5).                                            |
| 13             | 4096          | 000**1**000000000000 | Acceptor vector has a HWE-P value < 1e-4.                                             |
| 14             | 8192          | 00**1**0000000000000 | Donor vector has a HWE-P value < 1e-4.                                                |

## Viewing and manipulating output data
It is possible to filter `two` output data by: 
1) either start or end contig e.g. `chr1`, 
2) position in that contig e.g. `chr1:10e6-20e6`; 
3) have a particular contig mapping e.g. `chr1,chr2`; 
4) interval mapping in both contigs e.g. `chr1:10e3-10e6,chr2:0-10e6`

```bash
tomahawk view -i 1kgp3_chr6_1_45.two -I 6:10e3-10e6,6:0-10e6 -H | head -n 6
```


Viewing all data
```bash
tomahawk view -i 1kgp3_chr6_1_45.two -H | head -n 6
```

| FLAG | CHROM_A | POS_A  | CHROM_B | POS_B  | REF_REF | REF_ALT | ALT_REF | ALT_ALT | D           | DPrime   | R        | R2       | P           | ChiSqModel | ChiSqTable |
|------|---------|--------|---------|--------|---------|---------|---------|---------|-------------|----------|----------|----------|-------------|------------|------------|
| 2059 | 6       | 89572  | 6       | 214654 | 4999    | 6       | 0       | 3       | 0.000597965 | 1        | 0.577004 | 0.332934 | 4.01511e-09 | 1667.33    | 0          |
| 2059 | 6       | 89573  | 6       | 214654 | 4999    | 6       | 0       | 3       | 0.000597965 | 1        | 0.577004 | 0.332934 | 4.01511e-09 | 1667.33    | 0          |
| 2059 | 6       | 122855 | 6       | 214654 | 4988    | 17      | 0       | 3       | 0.000596649 | 1        | 0.38664  | 0.149491 | 5.44908e-08 | 748.648    | 0          |
| 3    | 6       | 143500 | 6       | 212570 | 4421    | 387     | 82      | 118     | 0.0195352   | 0.54402  | 0.331324 | 0.109775 | 1.53587e-69 | 549.755    | 0          |
| 3    | 6       | 143500 | 6       | 213499 | 4419    | 387     | 84      | 118     | 0.0194949   | 0.537523 | 0.329068 | 0.108286 | 7.57426e-69 | 542.296    | 0          |

Slicing a range
```bash
tomahawk view -i 1kgp3_chr6_1_45.two -H -I 6:5e6-6e6 | head -n 6
```

| FLAG | CHROM_A | POS_A   | CHROM_B | POS_B  | REF_REF | REF_ALT | ALT_REF | ALT_ALT | D           | DPrime   | R        | R2       | P            | ChiSqModel | ChiSqTable |
|------|---------|---------|---------|--------|---------|---------|---------|---------|-------------|----------|----------|----------|--------------|------------|------------|
| 2063 | 6       | 5022893 | 6       | 73938  | 5000    | 7       | 0       | 1       | 0.000199362 | 1        | 0.353306 | 0.124825 | 0.00159744   | 625.125    | 0          |
| 2063 | 6       | 5020564 | 6       | 89339  | 5000    | 7       | 0       | 1       | 0.000199362 | 1        | 0.353306 | 0.124825 | 0.00159744   | 625.125    | 0          |
| 7    | 6       | 5018459 | 6       | 156100 | 1150    | 1283    | 388     | 2187    | 0.0804323   | 0.509361 | 0.348864 | 0.121706 | 1.17649e-138 | 609.504    | 0          |
| 7    | 6       | 5018691 | 6       | 156100 | 861     | 811     | 677     | 2659    | 0.0693918   | 0.339199 | 0.31898  | 0.101748 | 1.10017e-109 | 509.554    | 0          |
| 7    | 6       | 5018910 | 6       | 156100 | 861     | 811     | 677     | 2659    | 0.0693918   | 0.339199 | 0.31898  | 0.101748 | 1.10017e-109 | 509.554    | 0          |

Slicing matches in B string
```bash
tomahawk view -i 1kgp3_chr6_1_45.two -H -I 6:5e6-6e6,6:5e6-6e6 | head -n 6
```

| FLAG | CHROM_A | POS_A   | CHROM_B | POS_B   | REF_REF | REF_ALT | ALT_REF | ALT_ALT | D           | DPrime   | R        | R2       | P           | ChiSqModel | ChiSqTable |
|------|---------|---------|---------|---------|---------|---------|---------|---------|-------------|----------|----------|----------|-------------|------------|------------|
| 3    | 6       | 5000012 | 6       | 5073477 | 4988    | 9       | 5       | 6       | 0.0011915   | 0.544089 | 0.465743 | 0.216917 | 1.05037e-13 | 1086.32    | 0          |
| 1027 | 6       | 5000160 | 6       | 5072168 | 5000    | 2       | 4       | 2       | 0.000398404 | 0.4994   | 0.407677 | 0.166201 | 7.1708e-06  | 832.333    | 0          |
| 1027 | 6       | 5000160 | 6       | 5078340 | 5000    | 2       | 4       | 2       | 0.000398404 | 0.4994   | 0.407677 | 0.166201 | 7.1708e-06  | 832.333    | 0          |
| 3    | 6       | 5000482 | 6       | 5079621 | 4993    | 6       | 2       | 7       | 0.0013931   | 0.777199 | 0.64641  | 0.417846 | 3.9492e-18  | 2092.57    | 0          |
| 3    | 6       | 5000482 | 6       | 5082227 | 4994    | 6       | 1       | 7       | 0.00139362  | 0.874675 | 0.685808 | 0.470333 | 8.78523e-19 | 2355.43    | 0          |

As alluded to in the [sorting section](#sorting-output-files), sorted files have 
generally much faster query times. We can demonstrate this difference by using the
sorted and unsorted data from the sliding window example above.
```bash hl_lines="3"
$ time tomahawk view -i 1kgp3_chr6_4mb.two -H -I 6:1e6-5e6,6:1e6-5e6 > /dev/null

real	3m1.029s
user	2m41.929s
sys	    0m5.682s
```

Sorted
```bash hl_lines="3"
$ time tomahawk view -i 1kgp3_chr6_4mb_sorted.two -H -I 6:1e6-5e6,6:1e6-5e6 > /dev/null

real	0m38.167s
user	0m37.579s
sys	    0m0.192s
```

Viewing `ld` data from the binary `two` file format and filtering out lines with a
Fisher's exact test P-value < 1e-4, minor haplotype frequency < 5 and have
both markers on the same contig (bit `2`)

```bash
tomahawk view -i file.two -P 1e-4 -a 5 -f 2
```

Example output

| FLAG | CHROM_A | POS_A   | CHROM_B | POS_B   | REF_REF | REF_ALT | ALT_REF | ALT_ALT | D             | Dprime | R          | R2         | P             | ChiSqModel | ChiSqTable |
|------|---------|---------|---------|---------|---------|---------|---------|---------|---------------|--------|------------|------------|---------------|------------|------------|
| 15   | 20      | 1314874 | 20      | 2000219 | 5002    | 5       | 0       | 1       | 0.00019944127 | 1      | 0.4080444  | 0.16650023 | 0.0011980831  | 0          | 833.83313  |
| 15   | 20      | 1315271 | 20      | 1992301 | 5005    | 2       | 0       | 1       | 0.00019956089 | 1      | 0.57723492 | 0.33320019 | 0.00059904153 | 0          | 1668.6665  |
| 15   | 20      | 1315527 | 20      | 1991024 | 5004    | 0       | 3       | 1       | 0.00019952102 | 1      | 0.49985018 | 0.2498502  | 0.00079872204 | 0          | 1251.2498  |
| 15   | 20      | 1315763 | 20      | 1982489 | 5006    | 0       | 1       | 1       | 0.00019960076 | 1      | 0.70703614 | 0.49990013 | 0.00039936102 | 0          | 2503.4999  |
| 15   | 20      | 1315807 | 20      | 1982446 | 5004    | 3       | 0       | 1       | 0.00019952102 | 1      | 0.49985018 | 0.2498502  | 0.00079872204 | 0          | 1251.2498  |

Example unphased output. Notice that estimated haplotype counts are now floating values and the bit-flag 1 is not set.

| FLAG | CHROM_A | POS_A   | CHROM_B | POS_B   | REF_REF   | REF_ALT   | ALT_REF       | ALT_ALT    | D             | Dprime     | R          | R2         | P             | ChiSqModel    | ChiSqTable |
|------|---------|---------|---------|---------|-----------|-----------|---------------|------------|---------------|------------|------------|------------|---------------|---------------|------------|
| 46   | 20      | 1882564 | 20      | 1306588 | 5003.9996 | 2.0003999 | 1.0003999     | 0.99960008 | 0.00019936141 | 0.49950022 | 0.40779945 | 0.1663004  | 0.0011978438  | 0.0011996795  | 832.83241  |
| 46   | 20      | 1895185 | 20      | 1306588 | 5004.9998 | 1.0001999 | 1.0001999     | 0.99980011 | 0.0001994811  | 0.49970025 | 0.49970022 | 0.24970032 | 0.00079864228 | 0.00069960005 | 1250.4992  |
| 46   | 20      | 1306588 | 20      | 1901581 | 5003.9996 | 2.0003999 | 1.0003999     | 0.99960008 | 0.00019936141 | 0.49950022 | 0.40779945 | 0.1663004  | 0.0011978438  | 0.0011996795  | 832.83241  |
| 46   | 20      | 1901581 | 20      | 1306588 | 5003.9996 | 2.0003999 | 1.0003999     | 0.99960008 | 0.00019936141 | 0.49950022 | 0.40779945 | 0.1663004  | 0.0011978438  | 0.0011996795  | 832.83241  |
| 46   | 20      | 1306649 | 20      | 1885268 | 5006      | 1         | 1.2656777e-08 | 0.99999999 | 0.00019960076 | 1          | 0.70703614 | 0.49990013 | 0.00039936102 | 0.00039969285 | 2503.4999  |

Example output for forced phased math in fast mode. Note that only the `REF_REF` count is available, the fast math bit-flag is set, and all P and Chi-squared CV values are 0.

| FLAG | CHROM_A | POS_A   | CHROM_B | POS_B   | REF_REF | REF_ALT | ALT_REF | ALT_ALT | D            | Dprime     | R          | R2         | P | ChiSqModel | ChiSqTable |
|------|---------|---------|---------|---------|---------|---------|---------|---------|--------------|------------|------------|------------|---|------------|------------|
| 67   | 20      | 1345922 | 20      | 1363176 | 4928    | 0       | 0       | 0       | 0.011788686  | 0.98334378 | 0.86251211 | 0.74392718 | 0 | 0          | 0          |
| 67   | 20      | 1345922 | 20      | 1367160 | 4933    | 0       | 0       | 0       | 0.011800847  | 0.98336071 | 0.89164203 | 0.79502547 | 0 | 0          | 0          |
| 67   | 20      | 1345958 | 20      | 1348347 | 4944    | 0       | 0       | 0       | 0.0092644105 | 0.97890127 | 0.85316211 | 0.72788554 | 0 | 0          | 0          |
| 67   | 20      | 1345958 | 20      | 1354524 | 4938    | 0       | 0       | 0       | 0.0092493389 | 0.86871892 | 0.80354655 | 0.64568704 | 0 | 0          | 0          |
| 75   | 20      | 1345958 | 20      | 1356626 | 4945    | 0       | 0       | 0       | 0.0033518653 | 0.99999994 | 0.51706308 | 0.26735422 | 0 | 0          | 0          |


## Aggregating and visualizing datasets

Tomahawk generally output many millions to many hundreds of millions to billions
of output linkage disequilibrium (LD) associations generated from many millions
of input SNVs. It is technically very challenging to visualize such large
datasets. Read more about [aggregation](aggregation.md) in Tomahawk.

| Aggregation | Description                                                       |
|-------------|-------------------------------------------------------------------|
| `R`           | Pearson correlation coefficient                                   |
| `R2`          | Squared pearson correlation coefficient                           |
| `D`           | Coefficient of linkage disequilibrium                             |
| `Dprime`      | Scaled coefficient of linkage disequilibrium                      |
| `Dp`          | Alias for `dprime`                                                |
| `P`           | Fisher's exact test P-value of the 2x2 haplotype contigency table |
| `Hets`        | Number of (0,1) or (1,0) associations                               |
| `Alts`        | Number of (1,1) associations                                |
| `Het`         | Alias for `hets`                                                  |
| `Alt`         | Alias for `alts`                                                  |

| Reduction | Description                             |
|-----------|-----------------------------------------|
| `Mean`      | Mean number of aggregate                |
| `Max`       | Largest number in aggregate bin         |
| `Min`       | Smallest number in aggregate bin        |
| `Count`     | Total number of records in a bin        |
| `N`         | Alias for `count`                       |
| `Total`     | Sum total of aggregated number in a bin |