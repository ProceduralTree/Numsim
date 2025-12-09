#!/usr/bin/env sh

for i in {1..4};do
    mpirun -n $i ../build/numsim_parallel ../settings.txt
   done 
