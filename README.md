![screenshot](tomahawk.png)
## Fast calculation of LD in large-scale cohorts
Synopsis

Marcus D. R. Klarqvist (<mk21@sanger.ac.uk>)

### Usage
### Installation instructions
Compiling Tomahawk with default parameters is straightforward.
```bash
git clone https://github.com/mklarqvist/Tomahawk_private
cd Tomahawk_private
cd build
make
```
By default, Tomahawk compiles using extremely aggressive optimization flags and
using native architecture instructions
(`-O3 -fno-strict-aliasing -march=native -mtune=native -ftree-vectorize -pipe
  -fomit-frame-pointer -flto -frename-registers -funroll-loops -fuse-linker-plugin`)
and internally compiles for the most recent SIMD-instruction set available.
This might result in additional effort when submitting array jobs to
computer farms/clouds with non-uniform hardware.

### Brief usage instructions
Tomahawk comprises five primary commands: `import`, `calc`, `view`, `sort`, and `index`
Executing `tomahawk <command>` gives a list of commands with brief descriptions

All primary Tomahawk commands operate on the binary Tomahawk `twk` and Totempole `twi` file
format. Interconversions between `twk` and `vcf`/`bcf` is supported through the
commands `import` for `vcf`/`bcf`->`twk` and `view` for `twk`->`vcf`/`bcf`. Linkage
disequilibrium data is written out in `two` and `toi` format. There is an option
to produce human-readable tab-delimited output for smaller datasets (`-N` flag in `calc`).

By default Tomahawk only operates on bi-allelic SNVs and as such filters out
indels and complex variants.
Importing a variant document (`vcf`/`bcf`) to Tomahawk requires the `import` command.
The following command line imports a `vcf` file and outputs `outPrefix.twk` and
`outPrefix.twk.twi`
```bash
tomahawk import -i file.vcf -o outPrefix
```

```bash
tomahawk calc -Bpdi file.twk -o - -a 5 -r 0.1 -P 0.1 -c 990 -C 1 -t 28 > output.two
```

 ```bash
 tomahawk view -i file.two
 ```

 ```bash
 tomahawk view -i file.twk
 ```

 ### License
 [MIT](LICENSE)
