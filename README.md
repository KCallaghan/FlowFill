# FlowFill

***FlowFill is a Fortran tool to fill or partially fill depressions in a landscape given a starting runoff value, in a way that conserves water mass.***

This algorithm is intended for pre-processing DEMs prior to running flow routing algorithms. Rather than indiscriminately filling all depressions, a starting runoff value is specified and depressions are fully or partially filled in a way that conserves water mass. Terrestrial water storage can also be assessed. 

When using FlowFill, please cite:

**Callaghan, K. L,. and A. D. Wickert (2019), [Computing water flow through complex landscapes, Part 1: Incorporating depressions in flow routing using FlowFill](https://www.earth-surf-dynam-discuss.net/esurf-2019-11/), *Earth Surf. Dynam. Discuss.*, doi:10.5194/esurf-2019-11.**

This code may only work on UNIX-based systems. GRASS GIS users may be interested in the r.flowfill extension, available at https://github.com/OSGeo/grass-addons/tree/master/grass7/raster/r.flowfill.

Please contact us if you have questions or suggestions!

## Required data inputs
The only required data file is topography in a .nc (NetCDF) format. The NetCDF file should have three variables: 'lat' for latitude, 'lon' for longitude, and 'value' for elevation.

Actual latitude and longitude are not required: define any rectangular grid. Note that this algorithm does not internally account for changing cell size in a geographic coordinate system. 

The filename, number of rows and columns in the file, runoff depth, and threshold value are specified at runtime. The threshold value is used to allow FlowFill to exit and save outputs once further computation will not make a significant difference to outputs. Appropriate values vary depending on the landscape used. We recommend running FlowFill with a very small threshold value and a small runoff amount and plotting the maximum amount of water moving per iteration (h_max) (FlowFill saves this to a text file) to aid in selection of the threshold. The threshold value represents the amount by which h_max should be allowed to vary from one iteration to the next. A distinctive plateau in h_max values is generally seen at the point where the threshold should be invoked.

## Dependencies

* The GNU Fortran compiler
* Open MPI
* NetCDF for Fortran

Install these on Ubuntu Linux using:
```
sudo apt install libopenmpi-dev
sudo apt install libnetcdff-dev
sudo apt install gfortran
```

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

Runoff and the topography file need to be specified at runtime, as well as the other user-selected variables. The easiest way to run FlowFill is to download the 'user_inputs.txt' file and change that values in that file as needed; download the 'run_me.sh' file and use this to run the code. 

The values that are specified in the 'user_inputs.txt' file are as follows:
1) Number of processors (integer value, 3 or larger). This will depend on how many processors you have available on your system and wish to dedicate to FlowFill. 
2) Name of compiled file to run, as specified in the compilation step above
3) Name of topography input file, in netcdf format
4) Columns in the topography input file
5) Rows in the topography input file
6) Threshold value - a small value that indicates when the code should stop running and equilibrium has been reached. The appropriate value varies on different input topographies. The threshold value indicates the amount of change in the maximum amount of water moving across the landscape from one iteration to the next. For example graphs on how this threshold visually appears, see the FlowFill paper. There is a distinct plateau in the maximum water moving that indicates what the appropriate threshold value shoule be. The threshold value is best selected by firstly choosing a very small value and running the code until the point where the amount of water moving from one iteration to the next has all but stopped. Use this information to select a more appropriate threshold for future runs of FlowFill on this topogrphy.
7) Output file prefix
8) Are you including a variable runoff file (Y/N) (as opposed to selecting a single contant value for runoff across the whole domain)
9) If N, Runoff (m)
10) If Y, name of variable runoff file 
11) Should ties be treated by selecting a preferential direction (PREF) or a random direction (RAND). Note that selecting RAND will reduce any biases introduced by ties, but will make the result non-deterministic (results may be different when performing multiple model runs). 
 
## Outputs

The program outputs two binary (.dat) files:
1. Filled landscape (i.e. original topography, plus depth of water in depressions)
2. Thickness of water (i.e. depth of water stored in depressions alone)

An additional output is a text file with runtimes and amounts of water still moving in the landscape every 500 iterations, along with other status messages. 
