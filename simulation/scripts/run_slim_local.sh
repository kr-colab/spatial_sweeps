#!/bin/bash


mkdir -p out/sim_$1\_s_$4\_Nw_$2
cd out/sim_$1\_s_$4\_Nw_$2
slim -d SI=$3 -d SM=$3 -d SD=$3 -d S=$4 -d TREETICK=$5 /home/crehmann/kernlab/spatial_sweeps/recipes/spatial_sweep.slim

#python ../../../slim_streamer.py --NW $3 --SC $4 --TT $5
