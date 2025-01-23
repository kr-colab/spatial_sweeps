#!/bin/bash

echo "file sim scoef Nw tick cutoff N_mutations frequency area distance outliers outlier_SNPdetected windowoutliers windowoutlier_SNPdetected windowoutlier_nonoutlier_SNPdetected" > power_analysis_sliding_cutoff.txt

i=1
#for f in $(find out -name '*_frequency_area.txt'); do
#    for cutoff in 0.01 0.05 0.1 0.2 0.3 0.4 0.5; do
while read f; do
        f="${f#< }"
        IFS='_'; read -ra FA <<< $f
        sim=${FA[1]}
        scoef=${FA[3]}
        Nw=${FA[5]}
        Nw=${Nw%/*}
        tick=${FA[6]}
        echo "$f" "$sim" "$scoef" "$Nw" "$tick"
        sbatch run_power_analysis.batch "$f" "$sim" "$scoef" "$Nw" "$tick" "power_analysis_sliding_cutoff.txt"
        if (( $i % 10000 == 0 )); then
            until (( "$(squeue -u crehmann | wc -l) < 10" )); do
                sleep 2m
            done
        fi
        i=$(( $i + 1 ))

done <frequency_area_files.txt

