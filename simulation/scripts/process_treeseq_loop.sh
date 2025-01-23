#!/bin/bash

# reprocess tree sequences that were zipped...
#for f in $(find out -name 'simulation_*.trees.gzip'); do
#    IFS='_'; read -ra FA <<< $f
#    sim=${FA[1]}
#    scoef=${FA[3]}
#    Nw=${FA[5]}
#    Nw=${Nw%/*}
#    tick=${FA[6]}
#    while (( $(squeue -u crehmann | wc -l) > 1500 )); do
#        sleep 1
#    done
#    echo "$f" "$sim" "$scoef" "$Nw" "$tick"
#    mv $f ${f%??}
#    gunzip $f        
#    sbatch process_treeseq.batch "$f" "$sim" "$scoef" "$Nw" 
#done
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
