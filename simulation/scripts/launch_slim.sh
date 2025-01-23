#!/bin/bash

declare -A sigmas
sigmas[10]=0.3989422804014327
sigmas[100]=1.2615662610100802
sigmas[1000]=3.989422804014327

#SCOEF=0.1
SCOEF=0.1
TREETICK=1

#sim_10_s_0.1_Nw_10

for Nw in 10 100 1000; do
    for replicate in {1..100}; do
        if [ ! -d out/sim_$replicate\_s_$SCOEF\_Nw_$Nw ]; then
            sbatch run_slim_local.batch $replicate $Nw ${sigmas[$Nw]} $SCOEF $TREETICK
           #echo $replicate $Nw ${sigmas[$Nw]}
#        else echo $replicate $Nw ${sigmas[$Nw]}
        fi
    done
done

#SCOEF=0.01
SCOEF=0.01
TREETICK=2
for Nw in 10 100 1000; do
    for replicate in {1..100}; do
        if [ ! -d out/sim_$replicate\_s_$SCOEF\_Nw_$Nw ]; then
            sbatch run_slim_local.batch $replicate $Nw ${sigmas[$Nw]} $SCOEF $TREETICK
#        else echo $replicate $Nw ${sigmas[$Nw]}
        fi
    done
done

#SCOEF=0.001
SCOEF=0.001
TREETICK=5
for Nw in 10 100 1000; do
    for replicate in {1..100}; do
        if [ ! -d out/sim_$replicate\_s_$SCOEF\_Nw_$Nw ]; then
            sbatch run_slim_local.batch $replicate $Nw ${sigmas[$Nw]} $SCOEF $TREETICK
#        else echo $replicate $Nw ${sigmas[$Nw]}
        fi
    done
done

#SCOEF=0.0001
SCOEF=0.0001
TREETICK=10

for Nw in 10 100 1000; do
    for replicate in {1..100}; do
        if [ ! -d out/sim_$replicate\_s_$SCOEF\_Nw_$Nw ]; then
            sbatch run_slim_local.batch $replicate $Nw ${sigmas[$Nw]} $SCOEF $TREETICK
#        else echo $replicate $Nw ${sigmas[$Nw]}
        fi
    done
done
