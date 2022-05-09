#!/bin/bash

 #LLVM 15 
 module use /sw/summit/modulefiles/ums/stf010/Core
 module load llvm/15.0.0-20220412 cuda/11.4.2

 #LLVM 14
 #module load llvm/14.0.0-latest cuda

 make COMPILER=llvm_nv TS=32
