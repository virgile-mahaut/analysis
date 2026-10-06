#!/bin/bash                                                                                                                                                                                               

export USER="$(id -u -n)"
export LOGNAME=${USER}
export HOME=/sphenix/u/${LOGNAME}/macros/detectors/sPHENIX/
export MYINSTALL=/sphenix/user/bkimelman/sPHENIX/install/

source /opt/sphenix/core/bin/sphenix_setup.sh -n
source /opt/sphenix/core/bin/setup_local.sh $MYINSTALL

DIR=/sphenix/tg/tg01/jets/bkimelman/VandyDSTs_wEEC_3D_unfolding_kinematics_Apr27_2026/

pids=()
for i in 12 20 30 40 50 60; do
    ../run_hadd_unfold.sh $i $DIR response kFull &
    pids+=($!)
    ../run_hadd_unfold.sh $i $DIR response kHalf &
    pids+=($!)
done
 
# Wait for all background jobs and check for failures
failed=0
for pid in "${pids[@]}"; do
    if ! wait "$pid"; then
        echo "ERROR: job $pid failed" >&2
        failed=1
    fi
done
if (( failed )); then
    echo "One or more hadd jobs failed — aborting" >&2
    exit 1
fi

#for i in 12 20 30 40 50 60; do
    #../run_hadd_unfold.sh $i $DIR pairWeights kFull
#    ../run_hadd_unfold.sh $i $DIR response kFull
#    ../run_hadd_unfold.sh $i $DIR response kHalf
#done

hadd -f $DIR/response-all-fullClosure.root $DIR/response_Jet*-fullClosure.root
rm $DIR/response_Jet*-fullClosure.root

hadd -f $DIR/response-all-halfClosure.root $DIR/response_Jet*-halfClosure.root
rm $DIR/response_Jet*-halfClosure.root

hadd -f -j 8 -n 100 $DIR/data_measured-all.root $DIR/Data/data_measured_*.root

echo "done with all hadds"
