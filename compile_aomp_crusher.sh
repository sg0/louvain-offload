#!/bin/bash

 #LLVM 15 

 module load craype-accel-amd-gfx90a
 module use /sw/crusher/ums/compilers/modulefiles/
 #module load llvm/15.0.0-20220425
 module load llvm

 make COMPILER=llvm_amd TS=32
