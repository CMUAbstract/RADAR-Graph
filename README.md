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

**NOTE**: For the `HUBDUP` and `RADAR` optimizations, the number of hubs to be 
duplicated should be selected based on the LLC Capacity (i.e. the `llcCap` 
variable should be selected based on the platform's available LLC capacity --
we found `llcCap` = 0.9 * _LLCSz_ offers the best performance)

## Compilation

The above versions have been tested using only OPENMP. To compile 
applications with OPENMP for parallelization, set the OPENMP environment
variable. 

Running `python compile_all.py` in the `src` directory will build all 
the applications with different optimizations

**NOTE**: Applications tested with g++-6.3.0 on debian stretch

## Usage Instructions

The graph input formats are the same as specified by Ligra ([link](https://github.com/jshun/ligra#input-format-for-ligra-applications-and-the-ligra-encoder)) and GAP([link](https://github.com/sbeamer/gapbs#graph-loading))

For symmetric input graphs, adding the `-s` flag will lead to a more efficient execution

## Citation

If this work was useful to you, please consider citing our work. 

bibtex entry

`@inproceedings{Balaji:2019:CDD:3307681.3326609,
 author = {Balaji, Vignesh and Lucia, Brandon},
 title = {Combining Data Duplication and Graph Reordering to Accelerate Parallel Graph Processing},
 booktitle = {Proceedings of the 28th International Symposium on High-Performance Parallel and Distributed Computing},
 series = {HPDC '19},
 year = {2019},
 isbn = {978-1-4503-6670-0},
 location = {Phoenix, AZ, USA},
 pages = {133--144},
 numpages = {12},
 url = {http://doi.acm.org/10.1145/3307681.3326609},
 doi = {10.1145/3307681.3326609},
 acmid = {3326609},
 publisher = {ACM},
 address = {New York, NY, USA},
 keywords = {atomics, data duplication, graph processing, graph reordering, locality, power-law},
}`


## Contact

For bugs or any other information, please contact:

`vigneshb_AT_andrew_DOT_cmu_DOT_edu`

