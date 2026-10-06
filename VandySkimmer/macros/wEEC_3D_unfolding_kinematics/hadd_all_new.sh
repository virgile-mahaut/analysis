#!/bin/bash

export USER="$(id -u -n)"
export LOGNAME=${USER}
export HOME=/sphenix/u/${LOGNAME}/macros/detectors/sPHENIX/
export MYINSTALL=/sphenix/user/bkimelman/sPHENIX/install/

source /opt/sphenix/core/bin/sphenix_setup.sh -n
source /opt/sphenix/core/bin/setup_local.sh $MYINSTALL

#DIR=/sphenix/tg/tg01/jets/bkimelman/VandyDSTs_wEEC_3D_unfolding_kinematics_Apr30_2026_0p25/
DIR=/sphenix/tg/tg01/jets/bkimelman/VandyDSTs_wEEC_3D_unfolding_kinematics_Apr30_2026_${1}/
run_pipeline() {
    local kMode=$1
    local pids=()
    local failed=0

    echo "Starting pipeline for $kMode"

    # --- step 1: launch jobs in parallel ---
    for i in 12 20 30 40 50 60; do
        ../run_hadd_unfold.sh $i $DIR response $kMode &
        pids+=($!)
    done

    # --- wait for step 1 jobs ---
    for pid in "${pids[@]}"; do
        if ! wait "$pid"; then
            echo "ERROR: $kMode job $pid failed" >&2
            failed=1
        fi
    done

    # --- only continue if successful ---
    if (( failed )); then
        echo "$kMode pipeline failed — skipping final hadd" >&2
        return 1
    fi

    echo "$kMode step1 complete, starting final hadd"

    # --- step 2: mode-specific hadd ---
    if [[ "$kMode" == "kFull" ]]; then
        hadd -f "$DIR/response-all-fullClosure.root" "$DIR"/response_Jet*-fullClosure.root
        rm "$DIR"/response_Jet*-fullClosure.root
    elif [[ "$kMode" == "kHalf" ]]; then
        hadd -f "$DIR/response-all-halfClosure.root" "$DIR"/response_Jet*-halfClosure.root
        rm "$DIR"/response_Jet*-halfClosure.root
    elif [[ "$kMode" == "kData" ]]; then
        hadd -f "$DIR/response-all-dataClosure.root" "$DIR"/response_Jet*-dataClosure.root
        rm "$DIR"/response_Jet*-dataClosure.root    
    elif [[ "$kMode" == "kVtx" ]]; then
        hadd -f "$DIR/response-all-vtx.root" "$DIR"/response_Jet*-vtx.root
        rm "$DIR"/response_Jet*-vtx.root
    fi

    echo "$kMode pipeline complete"
}

# --- launch both pipelines in parallel ---
run_pipeline kFull &
pid_full=$!

: '
run_pipeline kHalf &
pid_half=$!

run_pipeline kData &
pid_data=$!
'
# --- wait for both pipelines and track failures ---
fail_any=0

if ! wait "$pid_full"; then
    echo "kFull pipeline failed" >&2
    fail_any=1
fi

: '
if ! wait "$pid_half"; then
    echo "kHalf pipeline failed" >&2
    fail_any=1
fi

if ! wait "$pid_data"; then
    echo "kData pipeline failed" >&2
    fail_any=1
fi
'

# --- only run final hadd if both succeeded ---
if (( fail_any )); then
    echo "One or more pipelines failed — aborting final hadd" >&2
    exit 1
fi


echo "Running final data hadd"


#hadd -f -j 8 -n 100 "$DIR/data_measured-all.root" "$DIR"/Data/data_measured_*.root

echo "done with all hadds"