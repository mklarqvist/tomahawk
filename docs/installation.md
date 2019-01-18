# Installing Tomahawk

## Pre-requisites
* Unix-based operating system.
* Super-user priviledges if installing system-wide.

!!! warning "Requirement"
    
    Compiling and running Tomahawk requires a CPU with at least SSE4.2 available.
    This is a computational requirement stemming from the use of machine-optimized
    instructions called Single instruction, multiple data ([SIMD](https://en.wikipedia.org/wiki/SIMD)). This processing
    paradigm enables us to perform operations on multiple datapoints simultaneously.
    The ([Streaming SIMD Extensions](https://en.wikipedia.org/wiki/Streaming_SIMD_Extensions)) SSE4.2 instruction set, and later, describes 
    the instructions for performing these operations and are embedded in the CPU.

### Installation instructions
For modern x86-64 CPUs with `SSE4.2` or later, just type `make`. If you see
compilation errors, you most likely do not have `SSE4.2`. At the present time,
we do not support non-x86-64 CPUs or old CPU architecture.
```bash
git clone --recursive https://github.com/mklarqvist/tomahawk
cd tomahawk
make
```

### Debug mode
If you are extending upon Tomahawk or debugging, we provide a `DEBUG` flag to
build with all warnings triggered and with debug symbols enabled. See the main
[makefile](Makefile) for more information.
```
make DEBUG=true
```