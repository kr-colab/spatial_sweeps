#!/bin/bash

awk ' { print $1 } ' out/power_analysis.txt > finished_power_analysis.txt


diff <(sort frequency_area_files.txt) <(sort finished_power_analysis.txt) | grep "<" > unfinished_power_analysis.txt

