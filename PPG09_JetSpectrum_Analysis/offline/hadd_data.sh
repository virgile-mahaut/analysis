#!/bin/bash

radius_index=(2 3 4 5 6 8)

for ((r=0; r<${#radius_index[@]}; r++)); do
    rm output_data_r0${radius_index[$r]}.root
    hadd -j -k output_data_r0${radius_index[$r]}.root output_data/output_r0${radius_index[$r]}_*.root
done
