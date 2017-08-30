![screenshot](tomahawk.png)
## Fast calculation of LD in large-scale cohorts
### Synopsis

Marcus D. R. Klarqvist (<mk21@sanger.ac.uk>)

The current specification (v.0) is available [TWKv0](spec/TWKv0.pdf)

### Installation instructions
Compiling Tomahawk with default parameters is straightforward.
```bash
git clone https://github.com/mklarqvist/Tomahawk_private
cd Tomahawk_private
cd build
make
```
By default, Tomahawk compiles using extremely aggressive optimization flags and with
native architecture-specific instructions
(`-O3 -fno-strict-aliasing -march=native -mtune=native -ftree-vectorize -pipe
  -fomit-frame-pointer -flto -frename-registers -funroll-loops -fuse-linker-plugin`)
and internally compiles for the most recent SIMD-instruction set available.
This might result in additional effort when submitting jobs to
computer farms/clouds with a hardware architecture that is different from the compiled target.

### Brief usage instructions
Tomahawk comprises five primary commands: `import`, `calc`, `view`, `sort`, `index`,
`stats` and `concat`.
Executing `tomahawk <command>` gives a list of commands with brief descriptions

All primary Tomahawk commands operate on the binary Tomahawk `twk` and Totempole `twi` file
format. Interconversions between `twk` and `vcf`/`bcf` is supported through the
commands `import` for `vcf`/`bcf`->`twk` and `view` for `twk`->`vcf`. Linkage
disequilibrium data is written out in `two` and `toi` format. There is an option
to produce human-readable tab-delimited output for smaller datasets (`-N` flag in `calc`).

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
tomahawk calc -Bpdi file.twk -o - -a 5 -r 0.1 -P 0.1 -c 990 -C 1 -t 28 > output.two
```

### Converting between file formats and filtering
Viewing LD data from the binary `two` file format and filtering out lines with a
Fisher's exact test P-value < 1e-4, minor haplotype frequency < 5 and have
FLAG bits `4` set and map to any of these regions
 ```bash
 tomahawk view -i file.two -P 1e-4 -a 5 -f 4 '000001F|quiver:10e3-10e6,000004F|quiver:0-10e6' '000006F|quiver' '000007F|quiver:000009F|quiver'
 ```

 Viewing a `two` file
 ```bash
tomahawk view -Nhi file.two
 ```

Converting a `twk` file to `vcf`
 ```bash
 tomahawk view -HG -i file.twk -o file.vcf
 ```

 ### License
 [MIT](LICENSE)
