# Installing Tomahawk

# Pre-requisites
* Unix-based operating system.
* Super-user priviledges if installing system-wide.

!!! warning "Requirement"
    
    Compiling and running Tomahawk requires a CPU with at least SSE4.2 available.
    This is a computational requirement stemming from the use of machine-optimized
    instructions called Single instruction, multiple data (SIMD). This processing
    paradigm enables us to perform operations on multiple datapoints simultaneously.
    The (Streaming SIMD Extensions) SSE4.2 instruction set, and later, describes 
    the instructions for performing these operations and are embedded in the CPU.