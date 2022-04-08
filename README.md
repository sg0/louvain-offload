--------
Contents
--------
This folder contains codes for implementing Louvain method for 
graph community detection/clustering. 

This code requires an OpenMP runtime and a and C++11 compliant 
compiler for building. Experimental OpenMP offloading support is 
included (must pass -DUSE_OMP_OFFLOAD macro).

Please contact the following for any queries or support:
Sayan Ghosh, PNNL (sg0 at pnnl dot gov)

-------
Compile
-------
Just invoking `make should build the program, without any
changes made to the Makefile. Please pass -DUSE_OMP_OFFLOAD
and specific compiler options.

If you are running the code on a multi-socket system, pass
-DGRAPH_FT_LOAD=<x> while building where x == #sockets or 
#NUMA-nodes (it is 1 by default) to leverage first-touch access.

-----
Input
-----
We require graphs to be converted from its native format to a binary format.
The binary converter is part of another application, please follow the 
instructions for using Vite for file conversion: https://github.com/Exa-Graph/vite

Please note, this code requires weighted graphs (weight==1 is fine, outcome unknown 
for negative weights). There is an option in the Vite converter to ignore weights in 
graphs, which makes all edge weights zero. Such a graph won't be processed in
this case.

-------
Execute
-------
./matching_omp -f karate.bin

Apart from using external file, it is possible to generate
in-memory a random geometric graph in a distributed fashion.

Possible options (can be combined):

1. -f <bin-file>   : Specify input binary file after this argument. 
2. -n <vertices>   : Pass total number of vertices of the generated graph.
3. -l              : Use distributed LCG for randomly choosing edges. If this option 
                     is not used, we will use C++ random number generator (using 
                     std::default_random_engine).
4. -p <percent>    : Specify percent of overall edges to be randomly generated between
                     processes.
5. -w              : Use Euclidean distance as edge weight. If this option is not used,
                     edge weights are considered as 1.0. Generate edge weight uniformly 
                     between (0,1) if Euclidean distance is not available (applicable to 
                     randomly generated edges).  
6. -u              : Unit edge weights, i.e., make edge weights 1. This can potentially 
                     change outcome of the program, as original edge weights will be 
                     ignored when this option is used. If you are unsure about the weights
                     of the original graph and want to check a sample, invoke 
                     g->print_preview() on main (if nothing is shown, there is some issue 
                     with the binary graph itself). 
7. -h              : Print sample instances to run this code.                
