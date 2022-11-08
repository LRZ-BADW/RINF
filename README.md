# RINF Memory Bandwidth Benchmark 

## Introduction
This benchmark has been developed by Matthias Brehm and Reinhold Bader of LRZ.  

It tests the floating point and integer performance of various one-dimensional loop kernels. 
Depending on the access pattern, various aspects of the processor architecture and the memory hierarchy are tested. 
For example the `Double Precision Vector Triad`
	A(I) = B(I) * C(I) + D(I)
evaluated as 3 loads, 1 store and 2 floating-point operations is important (a potential load for “read for ownership” is not accounted).
The aggregate performance of a node is assessed. Scaling the benchmark for more than one node has the only purpose of assessing performance variations across cores or nodes.
Two subcases are considered:
•	rinf1.small.conf: probably best suited for filling the nodes only with (many) MPI tasks, without multithreading
•	rinf1.large.conf: probably best suited for filling the nodes with just one
In measures the memory bandwith within the nodes of a cluster. 

The benchmark code is written in standard Fortran, with some infrastructure as well as kernel code in standard C. 
Support for C interoperability and dynamic memory management from the Fortran 2003 standard is required. 
The MPI variant of the code presumes the availability of a standard-conforming MPI implementation 
including the Fortran module mpi; but only functionality from MPI-1 is needed. 

## License
This code is provided under the GNU General Public License v3.0.
Please find its terms in the accompanying file LICENSE.

## Build
To build the code, run 

    cd src; make

Resulting binaries are:
* rinf1.small.exe
* rinf1.large.exe

## Run

See src/JOB as example of job script

Example reference outputs are provided in the subdirectory REFERENCE-OUTPUT.

