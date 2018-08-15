# GP-WENO


Some Notes:
* Parallelization is done using Fortran Coarrays. Works with GNU compiler only if run as serial or with the OpenCoarrays library (http://www.opencoarrays.org/). Works with intel compiler using the makefile.intel
* HDF5 is needed for I/O
* Runtime Parameters can be set in the slug.init file
