#!/bin/bash


#find out -mindepth 1 -maxdepth 1 -type d '!' -exec bash -c '[[ -n $(find "{}" -iname "simulation_*_frequency_area.txt") ]]' \; -print > redo_sims.txt

#find out -mindepth 1 -maxdepth 1 -type d '!' -exec bash -c '[[ -n $(find "{}" -iname "simulation_*.trees*") ]]' \; -print > redo_sims.txt


find out -mindepth 1 -maxdepth 1 -type d '!' -exec bash -c '[[ -n $(find "{}" -iname "simulation.trees") ]]' \; -print > redo_sims.txt

while read d; do
    rm $d\/*
    rmdir $d
done <redo_sims.txt

   
