#!/bin/bash
# ─────────────────────────────────────────────────────────────────────────────
#  run_makeVtxWeights.sh  —  runs drawClosure.C
#
#  Usage:  ./run_makeVtxWeights.sh
#    DIR  — output directory
# ─────────────────────────────────────────────────────────────────────────────

export USER="$(id -u -n)"
export LOGNAME=${USER}
export HOME=/sphenix/u/${LOGNAME}/macros/detectors/sPHENIX/
export MYINSTALL=/sphenix/user/bkimelman/sPHENIX/install/

source /opt/sphenix/core/bin/sphenix_setup.sh -n
source /opt/sphenix/core/bin/setup_local.sh $MYINSTALL

DIR=/sphenix/tg/tg01/jets/bkimelman/VandyDSTs_wEEC_3D_unfolding_kinematics_Apr27_2026/

root -b -q "makeVtxWeights.C(\"${DIR}\")"

echo "done"