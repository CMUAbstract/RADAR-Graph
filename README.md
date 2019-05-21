# RADAR-Graph

Repository associated with _Combining Data Duplication and Graph Reordering to 
Accelerate Parallel Graph Processing_ **HPDC19**

## Optimizations included

The following optimizations (along with baseline push-only versions) can 
be found in `src`:
* Push-Pull
* HUBDUP
* DegreeSorting
* RADAR (_combination of the above two optimizations_)

The above codes were used for generating the main plots in the paper. Source
code for other experiments will be added soon...

## Compilation

The above versions have been tested using only OPENMP. To compile 
applications with OPENMP for parallelization, set the OPENMP environment
variable. 

Running `python compile_all.py` in the `src` directory will build all 
the applications with different optimizations

**NOTE**: Applications tested with g++-6.3.0 on debian stretch

## Contact

For bugs or any other information, please contact:

`vigneshb_AT_andrew_DOT_cmu_DOT_edu`

