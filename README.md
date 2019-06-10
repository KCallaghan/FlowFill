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
