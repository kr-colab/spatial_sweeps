#!/bin/bash

i=1
for f in $(find out/ -name 'simulation_*.trees'); do
    #fa=${f/.trees/_frequency_area.txt}
    #if [ ! -f "$fa" ]; then
    IFS='_'; read -ra FA <<< $f
    sim=${FA[1]}
    scoef=${FA[3]}
    Nw=${FA[5]}
    Nw=${Nw%/*}
    tick=${FA[6]}
    echo "$f" "$sim" "$scoef" "$Nw" "$tick"
    sbatch process_treeseq.batch "$f" "$sim" "$scoef" "$Nw" 
    if (( $i % 2000 == 0 )); then
        sleep 5m
    fi
    i=$(( $i + 1 ))
    echo "$i"
done
