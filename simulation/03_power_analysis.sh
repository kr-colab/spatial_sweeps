#!/bin/bash

echo "file sim scoef Nw tick N_mutations frequency area distance outliers outlier_SNPdetected windowoutliers windowoutlier_SNPdetected windowoutlier_nonoutlier_SNPdetected" > out/power_analysis.txt

i=1
for f in $(find out -name '*_frequency_area.txt'); do
    IFS='_'; read -ra FA <<< $f
    sim=${FA[1]}
    scoef=${FA[3]}
    Nw=${FA[5]}
    Nw=${Nw%/*}
    tick=${FA[6]}
    echo "$f" "$sim" "$scoef" "$Nw" "$tick"
    sbatch run_power_analysis.batch "$f" "$sim" "$scoef" "$Nw" "$tick" "out/power_analysis.txt"
    if (( $i % 10000 == 0 )); then
        sleep 5m
    fi
    i=$(( $i + 1 ))
done 

