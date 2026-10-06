#!/bin/bash

# Check if exactly three arguments are passed
if [ "$#" -ne 3 ]; then
    echo "Error: Missing required parameters."
    echo "Usage: $0 <do_reweighting (0 or 1)> <mrad_option (0 or 1)> <interaction_index (1 or 2)>"
    exit 1
fi

do_reweighting=$1
mrad_option=$2
interaction_index=$3

# Check if do_reweighting is strictly 0 or 1
if [ "$do_reweighting" != "0" ] && [ "$do_reweighting" != "1" ]; then
    echo "Error: do_reweighting parameter must be exactly 0 or 1."
    echo "Usage: $0 <do_reweighting (0 or 1)> <mrad_option (0 or 1)> <interaction_index (1 or 2)>"
    exit 1
fi

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

mkdir -p condor_sim/condor_sim_${suffix}
mkdir -p output_sim/output_sim_${suffix}

rm -f condor_sim/condor_sim_${suffix}/*
rm -f output_sim/output_sim_${suffix}/*

runindex=28
segperjob=50
radius=(0.4)
runtypes=("Jet5GeV" "Jet12GeV" "Jet20GeV" "Jet30GeV" "Jet40GeV" "Jet50GeV")

csv_file="ListCondorRunSim_${suffix}.csv"
sub_file="RunCondorSim.sub"

> "$csv_file"

# Set treepath based on interaction_index
if [ "$interaction_index" == "1" ]; then
    treepath="/sphenix/user/hanpuj/Dataset/Run${runindex}"
else
    treepath="/sphenix/user/hanpuj/Dataset/Run${runindex}_double"
fi

for runtype in "${runtypes[@]}"; do
    
    totalseg=$(ls "$treepath"/output_${runtype}_*.root 2>/dev/null | wc -l)
    total_condorjobs=$(( (totalseg + segperjob - 1) / segperjob ))
    
    echo "Processing Run type: ${runtype} (Total segments: ${totalseg})"

    for ((i=0; i<total_condorjobs; i++)); do
        if [[ $(((i+1)*segperjob)) -gt $totalseg ]]; then
            nseg=$((totalseg - i*segperjob))
        else
            nseg=$segperjob
        fi
        iseg=$((i*segperjob))

        for radius_value in "${radius[@]}"; do
            echo "$runtype,$nseg,$iseg,$radius_value,$do_reweighting,$mrad_option,$interaction_index,$mrad_str,$interaction_str" >> "$csv_file"
        done
    done
done

echo "Argument list formed: $csv_file"

condor_submit "$sub_file" list_file="$csv_file"

echo "All jobs submitted!"
