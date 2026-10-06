#!/bin/bash

source /sphenix/u/hj2552/setup.sh

Run="$1"
radius_value="$2"

echo "Run: ${Run}, radius: ${radius_value}"

root -l -q -b "analysis_data.C(${Run}, ${radius_value})"

echo "JOB COMPLETE!"
