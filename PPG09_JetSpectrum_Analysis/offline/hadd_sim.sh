#!/bin/bash

# Check if exactly two arguments are passed
if [ "$#" -ne 2 ]; then
    echo "Error: Missing required parameters."
    echo "Usage: $0 <mrad_option (0 or 1)> <interaction_index (1 or 2)>"
    exit 1
fi

mrad_option=$1
interaction_index=$2

# Check if mrad_option is strictly 0 or 1
if [ "$mrad_option" != "0" ] && [ "$mrad_option" != "1" ]; then
    echo "Error: mrad_option parameter must be exactly 0 or 1."
    echo "Usage: $0 <do_reweighting (0 or 1)> <mrad_option (0 or 1)> <interaction_index (1 or 2)>"
    exit 1
fi

# Check if interaction_index is strictly 1 or 2
if [ "$interaction_index" != "1" ] && [ "$interaction_index" != "2" ]; then
    echo "Error: interaction_index parameter must be exactly 1 or 2."
    echo "Usage: $0 <do_reweighting (0 or 1)> <mrad_option (0 or 1)> <interaction_index (1 or 2)>"
    exit 1
fi

# Set mrad string based on mrad_option
if [ "$mrad_option" == "0" ]; then
    mrad_str="0mrad"
else
    mrad_str="1p5mrad"
fi

# Set interaction string based on interaction_index
if [ "$interaction_index" == "1" ]; then
    interaction_str="single"
else
    interaction_str="double"
fi

# Construct suffix for directory and file names
suffix="${mrad_str}_${interaction_str}"

mkdir -p output_sim_hadd

rm output_sim_hadd/output_sim_${suffix}_*

radius_index=(4)

for ((r=0; r<${#radius_index[@]}; r++)); do
    hadd -j -k output_sim_hadd/output_sim_${suffix}_r0${radius_index[$r]}_MB.root output_sim/output_sim_${suffix}/output_r0${radius_index[$r]}_MB_*
    hadd -j -k output_sim_hadd/output_sim_${suffix}_r0${radius_index[$r]}_Jet5GeV.root output_sim/output_sim_${suffix}/output_r0${radius_index[$r]}_Jet5GeV_*
    hadd -j -k output_sim_hadd/output_sim_${suffix}_r0${radius_index[$r]}_Jet12GeV.root output_sim/output_sim_${suffix}/output_r0${radius_index[$r]}_Jet12GeV_*
    hadd -j -k output_sim_hadd/output_sim_${suffix}_r0${radius_index[$r]}_Jet20GeV.root output_sim/output_sim_${suffix}/output_r0${radius_index[$r]}_Jet20GeV_*
    hadd -j -k output_sim_hadd/output_sim_${suffix}_r0${radius_index[$r]}_Jet30GeV.root output_sim/output_sim_${suffix}/output_r0${radius_index[$r]}_Jet30GeV_*
    hadd -j -k output_sim_hadd/output_sim_${suffix}_r0${radius_index[$r]}_Jet40GeV.root output_sim/output_sim_${suffix}/output_r0${radius_index[$r]}_Jet40GeV_*
    hadd -j -k output_sim_hadd/output_sim_${suffix}_r0${radius_index[$r]}_Jet50GeV.root output_sim/output_sim_${suffix}/output_r0${radius_index[$r]}_Jet50GeV_*
    hadd -j -k output_sim_hadd/output_sim_${suffix}_r0${radius_index[$r]}_Jet60GeV.root output_sim/output_sim_${suffix}/output_r0${radius_index[$r]}_Jet60GeV_*
done
