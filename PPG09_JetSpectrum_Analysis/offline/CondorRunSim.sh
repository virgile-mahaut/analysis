#!/bin/bash

source /sphenix/u/hj2552/setup.sh

runtype="$1"
nSeg="$2"
iSeg="$3"
radius_value="$4"
do_reweighting="$5"
mrad_option="$6"
interaction_index="$7"

echo "Run type: ${runtype}, nSeg: ${nSeg}, iSeg: ${iSeg}, radius: ${radius_value}, do_reweighting: ${do_reweighting}, mrad_option: ${mrad_option}, interaction_index: ${interaction_index}"

root -l -q -b "analysis_sim.C(\"${runtype}\", ${nSeg}, ${iSeg}, ${radius_value}, ${do_reweighting}, ${mrad_option}, ${interaction_index})"

echo "JOB COMPLETE!"
