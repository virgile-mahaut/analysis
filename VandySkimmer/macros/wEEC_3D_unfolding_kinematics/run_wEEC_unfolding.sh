#!/bin/bash
# ─────────────────────────────────────────────────────────────────────────────
#  run_wEEC_unfolding.sh  —  runs wEEC_doUnfolding.C
#
#  Usage:  ./run_wEEC_unfolding.sh <DIR> <MODE> [NITER_WEEC]
#    DIR        — output directory
#    MODE       — kFull | kHalf | kData
#    NITER_WEEC — Bayes iterations for wEEC unfolding (default: 4)
#
#  Can be run independently of run_unfolding.sh — does not need jet_pT file.
#  For kData: uses fullClosure response for MC corrections,
#             data_measured-all.root for measured wEEC.
# ─────────────────────────────────────────────────────────────────────────────

export USER="$(id -u -n)"
export LOGNAME=${USER}
export HOME=/sphenix/u/${LOGNAME}/macros/detectors/sPHENIX/
export MYINSTALL=/sphenix/user/bkimelman/sPHENIX/install/

source /opt/sphenix/core/bin/sphenix_setup.sh -n
source /opt/sphenix/core/bin/setup_local.sh $MYINSTALL

DIR=/sphenix/tg/tg01/jets/bkimelman/VandyDSTs_wEEC_3D_unfolding_kinematics_Apr29_2026/
MODE=$1
NITER=${2:-4}

if [ -z "$DIR" ] || [ -z "$MODE" ]; then
    echo "Usage: $0 <DIR> <MODE> [NITER_WEEC]"
    exit 1
fi

case $MODE in
    kFull) LABEL="fullClosure" ;;
    kHalf) LABEL="halfClosure" ;;
    kData) LABEL="dataClosure" ;;
    *)
        echo "Unknown mode: $MODE  (expected kFull, kHalf, or kData)"
        exit 1
        ;;
esac

OUT_FILE="${DIR}/wEEC-${LABEL}.root"

if [ "$MODE" = "kData" ]; then
    RESP_FILE="${DIR}/response-all-${LABEL}.root"
    MEAS_FILE="${DIR}/data_measured-all.root"
    echo "Running wEEC_doUnfolding: mode=${MODE}, resp=${RESP_FILE}, meas=${MEAS_FILE}, out=${OUT_FILE}"
    root -b -q "wEEC_doUnfolding.C(\"${RESP_FILE}\",\"${OUT_FILE}\",${NITER},Mode::${MODE},\"${MEAS_FILE}\")"
else
    RESP_FILE="${DIR}/response-all-${LABEL}.root"
    echo "Running wEEC_doUnfolding: mode=${MODE}, resp=${RESP_FILE}, out=${OUT_FILE}"
    root -b -q "wEEC_doUnfolding.C(\"${RESP_FILE}\",\"${OUT_FILE}\",${NITER},Mode::${MODE})"
fi

echo "Done: ${OUT_FILE}"
