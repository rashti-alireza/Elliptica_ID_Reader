# Elliptica Initial Data Reader

This repository is a stand-alone library and provides an interface for an evolution code to
read the initial data of *Elliptica*.

Requirements:
- A C programming compiler, for instance, `gcc`.
- A compatible `openmp` library with the compiler.

Usage:
1. We choose the compiler at the `GNUmakefile` file. For example, a `gcc` compiler is 
chosen by `CC = gcc` and for the `openmp` library we set `CFLAGS += -fopenmp`.

2. `$ make`

3. To link this library to *an evolution code*, one needs to pass the path `-L/path/to/Elliptica_ID_Reader/lib` 
to the compiler during the compilation of the evolution code. 
Similarly, the compiler requires the path `-I/path/to/Elliptica_ID_Reader/include` 
to find the pertinent header file.

