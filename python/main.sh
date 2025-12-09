#!/usr/bin/env sh

for i in {1..8};do
    mpirun -n $i ../build/numsim_parallel ../settings.txt
   done 
