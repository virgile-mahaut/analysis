#!/bin/bash                                                                                                                                                                                               

export USER="$(id -u -n)"
export LOGNAME=${USER}
export HOME=/sphenix/u/${LOGNAME}/macros/detectors/sPHENIX/
export MYINSTALL=/sphenix/user/bkimelman/sPHENIX/install/

source /opt/sphenix/core/bin/sphenix_setup.sh -n
source /opt/sphenix/core/bin/setup_local.sh $MYINSTALL

jetSample=$1
segStart=$2
segStart=$((10#$segStart))
outDir=$3
mode=$4

root -b -q fillResponse.C\($jetSample,$segStart,\"$outDir\",Mode::$mode\)

echo "done with fillResponse job"
