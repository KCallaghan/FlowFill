# FlowFill

***FlowFill is a Fortran tool to fill or partially fill depressions in a landscape given a starting runoff value, in a way that conserves water mass.***

This algorithm is intended for pre-processing DEMs prior to running flow routing algorithms. Rather than indiscriminately filling all depressions, a starting runoff value is specified and depressions are fully or partially filled in a way that conserves water mass. Terrestrial water storage can also be assessed. 

This code may only work on UNIX-based systems. 

Please contact us if you have questions or suggestions! 

## Required data inputs
The only required data file is topography in a .nc (NetCDF) format. 

## Compilation

This code uses MPI and NetCDF libraries. Compile and run using MPI. 
An example of compilation is:

```
mpif90 -O3 surface_water.f90 -o Your_output -I/usr/include -L/usr/lib/x86_64-linux-gnu/ -lnetcdff -ffast-math -march=native -mtune=native
```
To use this compilation code on your computer, change the following to the appropriate paths:

* `-I/usr/include`
* `-L/usr/lib/x86_64-linux-gnu/`

## Running the code

This code will run in parallel and requires a minimum of 3 processors.
The runoff to be applied to the landscape should be specified at runtime using flag --runoff. 
An example of how to run this code with 1 m runoff and using 4 processors is:

```
 mpirun -np 4 ./Sangamon_test --runoff 1
 ```
 
## Outputs

The program outputs two binary (.dat) files:
1. Filled landscape (i.e. original topography, plus water in depressions)
2. Thickness of water (i.e. water stored in depressions alone)
