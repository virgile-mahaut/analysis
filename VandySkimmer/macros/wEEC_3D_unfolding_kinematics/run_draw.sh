#!/bin/bash
# ─────────────────────────────────────────────────────────────────────────────
#  run_draw.sh  —  runs drawClosure.C
#
#  Usage:  ./run_draw.sh <DIR> <MODE>
#    DIR  — output directory
#    MODE — kFull | kHalf | kData
#
#  Must be run after run_wEEC_unfolding.sh.
#  For closure modes the merged response file is passed so drawClosure can
#  load truth-tower wEEC histograms for the closure reference.
#  Creates DIR/Plots/ if it does not exist.
# ─────────────────────────────────────────────────────────────────────────────

export USER="$(id -u -n)"
export LOGNAME=${USER}
export HOME=/sphenix/u/${LOGNAME}/macros/detectors/sPHENIX/
export MYINSTALL=/sphenix/user/bkimelman/sPHENIX/install/

source /opt/sphenix/core/bin/sphenix_setup.sh -n
source /opt/sphenix/core/bin/setup_local.sh $MYINSTALL

#DIR=/sphenix/tg/tg01/jets/bkimelman/VandyDSTs_wEEC_3D_unfolding_kinematics_Apr28_2026/
DIR=/sphenix/tg/tg01/jets/bkimelman/VandyDSTs_wEEC_3D_unfolding_kinematics_Apr29_2026/
MODE=$1

if [ -z "$DIR" ] || [ -z "$MODE" ]; then
    echo "Usage: $0 <DIR> <MODE>"
    exit 1
fi

case $MODE in
    kFull) LABEL="fullClosure" ;;
    kHalf) LABEL="halfClosure" ;;
    kData) LABEL="dataClosure" ;;
    *)
        echo "Unknown mode: $MODE  (expected kFull, kHalf, or kData)"
        exit 1 ;;
esac

mkdir -p "${DIR}/Plots"

if [ "$MODE" = "kData" ]; then
    # kData: no truth reference needed, pass null for respFile
    RESP_FILE="${DIR}/response-all-${LABEL}.root"
    echo "Running drawClosure: mode=${MODE}, dir=${DIR}, resp=${RESP_FILE}"
    root -b -q "draw_misses_fakes.C(\"${DIR}\",\"${RESP_FILE}\",Mode::${MODE})"
    root -b -q "drawClosure.C(\"${DIR}\",\"${RESP_FILE}\",Mode::${MODE})"
else
    RESP_FILE="${DIR}/response-all-${LABEL}.root"
    if [ "$MODE" = "kFull" ]; then
        root -b -q "draw_3DClosure.C(\"${DIR}\")"
        root -b -q "drawResponse.C(\"${DIR}\",\"${RESP_FILE}\",Mode::${MODE},7,2,2,5)"
    fi
    root -b -q "draw_misses_fakes.C(\"${DIR}\",\"${RESP_FILE}\",Mode::${MODE})"
    echo "Running drawClosure: mode=${MODE}, dir=${DIR}, resp=${RESP_FILE}"
    root -b -q "drawClosure.C(\"${DIR}\",\"${RESP_FILE}\",Mode::${MODE})"
fi

echo "Done. Plots written to ${DIR}/Plots/"
