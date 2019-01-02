# FlowFill

***FlowFill is a Fortran tool to fill or partially fill depressions in a landscape given a starting runoff value, in a way that conserves water mass.***

This algorithm is intended for pre-processing DEMs prior to running flow routing algorithms. Rather than indiscriminately filling all depressions, a starting runoff value is specified and depressions are fully or partially filled in a way that conserves water mass. Terrestrial water storage can also be assessed. 

This code may only work on UNIX-based systems. 

Please contact us if you have questions or suggestions! 

## Required data inputs
The only required data file is topography in a .nc (NetCDF) format. The filename is specified at runtime. The NetCDF file should have three variables: 'lat' for latitude, 'lon' for longitude, and 'value' for elevation.

Actual latitude and longitude are not required: define any rectangular grid. Note that this algorithm does not internally account for changing cell size in a geographic coordinate system. 


## Compilation

This code uses MPI and NetCDF libraries. Compile and run using MPI. 
An example of compilation is:

```
mpif90 -O3 FlowFill.f90 -o Your_compiled_code -I/usr/include -L/usr/lib/x86_64-linux-gnu/ -lnetcdff -ffast-math -march=native -mtune=native
```
To use this compilation code on your computer, change the following to the appropriate paths:

* `-I/usr/include`
* `-L/usr/lib/x86_64-linux-gnu/`

## Running the code

This code will run in parallel and requires a minimum of 3 processors.

Runoff and the topography file need to be specified at runtime. Note that these must be specified in the correct order, first runoff and then the topography file.

An example of how to run this code with 1 m runoff, an input topography file called topo.nc, and using 4 processors is:

```
 mpirun -np 4 ./Your_compiled_code  1 topo.nc
 ```
 
## Outputs

The program outputs two binary (.dat) files:
1. Filled landscape (i.e. original topography, plus depth of water in depressions)
2. Thickness of water (i.e. depth of water stored in depressions alone)

An additional output is a text file with runtimes and amounts of water still moving in the landscape every 500 iterations, along with other status messages. 
