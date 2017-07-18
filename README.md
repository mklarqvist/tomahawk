![screenshot](tomahawk.png)
## Fast calculation of LD in large-scale cohorts


### Usage
### Building
```bash
cd build
make
```

### Examples
Tomahawk comprises five functions: `import`, `calc`, `view`, `sort`, and `index`

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
