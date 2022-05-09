#!/bin/bash


 module load craype-accel-amd-gfx90a
 module load cce rocm

 make COMPILER=cce TS=32
