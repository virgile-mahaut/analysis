#!/bin/bash
mkdir -p condor_data
mkdir -p output_data

rm condor_data/*
rm output_data/*

radius=(0.4)

> ListCondorRunData.csv

echo "Starting to form arguement list..."
for ((r=0; r<${#radius[@]}; r++)); do
    export radius_value=${radius[$r]}

    for run in $(cat listrunnumber_1p5mrad.txt); do
        echo "Run: ${run}"
        echo "$run,$radius_value" >> ListCondorRunData.csv
    done
done
echo "Arguement list formed!"

condor_submit RunCondorData.sub

echo "All jobs submitted!"
