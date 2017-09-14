![screenshot](tomahawk.png)
## Fast calculation of LD in large-scale cohorts
### Synopsis

Marcus D. R. Klarqvist (<mk21@sanger.ac.uk>)

The current specification (v.0) is available [TWKv0](spec/TWKv0.pdf)

### Installation instructions
Compiling Tomahawk with default parameters is straightforward.
```bash
git clone --recursive https://github.com/mklarqvist/Tomahawk
cd Tomahawk
cd build
make
```
By default, Tomahawk compiles using extremely aggressive optimization flags and with
native architecture-specific instructions
(`-march=native -mtune=native -ftree-vectorize -pipe -frename-registers -funroll-loops`)
and internally compiles for the most recent SIMD-instruction set available.
This might result in additional effort when submitting jobs to
computer farms/clouds with a hardware architecture that is different from the compiled target.

### Brief usage instructions
Tomahawk comprises five primary commands: `import`, `calc`, `view`, `sort`, and `concat`.
The functions `index` and `stats` are disabled at the moment.
Executing `tomahawk` gives a list of commands with brief descriptions and `tomahawk <command>`
gives detailed details for that command.

All primary Tomahawk commands operate on the binary Tomahawk `twk` and Totempole `twi` file
format. Interconversions between `twk` and `vcf`/`bcf` is supported through the
commands `import` for `vcf`/`bcf`->`twk` and `view` for `twk`->`vcf`. Linkage
disequilibrium data is written out in `two` and `toi` format.

### Importing to Tomahawk
By design Tomahawk only operates on bi-allelic SNVs and as such filters out
indels and complex variants. Tomahawk does not support mixed phasing of genotypes
in the same variant (e.g. `0|0`, `0/1`). If mixed phasing is found in a line,
all genotypes in that line are converted to unphased. Importing a variant document (`vcf`/`bcf`)
to Tomahawk requires the `import` command.
The following command line imports a `vcf` file and outputs `outPrefix.twk` and
`outPrefix.twk.twi` and filters out variants with >20% missingness and deviate
from Hardy-Weinberg equilibrium with a probability < 0.001
```bash
tomahawk import -i file.vcf -o outPrefix -m 0.2 -H 1e-3
```

### Import-extend
If you have split up your `vcf`/`bcf` files into multiple disjoint files (such as one per chromosome) it is possible to iteratively import and extend a `twk` file:
```bash
tomahawk import -i file.bcf -e extend.twk -m 0.2 -H 1e-3
```

### Calculating linkage disequilibrium
```bash
tomahawk calc -pdi file.twk -o output_prefix -a 5 -r 0.1 -P 0.1 -c 990 -C 1 -t 28
```

### Converting between file formats and filtering
Viewing LD data from the binary `two` file format and filtering out lines with a
Fisher's exact test P-value < 1e-4, minor haplotype frequency < 5 and have
FLAG bits `4` set
```bash
tomahawk view -i file.two -P 1e-4 -a 5 -f 4
 ```

It is possible to filter `two` output data by: 1) either start or end contig e.g.
`chr1`, 2) position in that contig `chr1:10e6-20e6`, or 3) or have a particular
contig mapping `chr1,chr2`, or 4) a particular regional mapping in both contigs
`chr1:10e3-10e6,chr2:0-10e6`
```bash
tomahawk view -i file.two `chr1:10e3-10e6,chr2:0-10e6`
 ```

Converting a `twk` file to `vcf`
 ```bash
tomahawk view -i file.twk -o file.vcf
```

### Sort `TWO` file
Partially sort `two` file in 500 MB chunks
```bash
tomahawk sort -i file.two -o partial.two -L 500
```

Perform k-way merge of partially sorted blocks
```bash
tomahawk sort -i partial.two -o sorted.two -M
```

### License
[MIT](LICENSE)
