# Tomahawk
## Fast calculation of LD in large-scale cohorts

![screenshot](tomahawk.png)

### Building
```bash
cd build
make
```

### Example
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
