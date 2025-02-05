#!/bin/bash

mkdir -p out/sim_$1\_s_$4\_Nw_$2
cd out/sim_$1\_s_$4\_Nw_$2
python ../../../slim_streamer.py --NW $3 --SC $4 --TT $5
