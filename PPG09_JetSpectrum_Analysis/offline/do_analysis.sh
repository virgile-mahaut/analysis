#!/usr/bin/env bash
# Run steps:
# sh run_condor_sim.sh 0 0 1; sh run_condor_sim.sh 0 0 2; sh run_condor_sim.sh 0 1 1; sh run_condor_sim.sh 0 1 2; sh run_condor_data.sh;
# sh hadd_sim.sh 0 1; sh hadd_sim.sh 0 2; sh hadd_sim.sh 1 1; sh hadd_sim.sh 1 2; sh hadd_data.sh;
# sh do_analysis.sh 0
# sh run_condor_sim.sh 1 0 1; sh run_condor_sim.sh 1 0 2; sh run_condor_sim.sh 1 1 1; sh run_condor_sim.sh 1 1 2;
# sh hadd_sim.sh 0 1; sh hadd_sim.sh 0 2; sh hadd_sim.sh 1 1; sh hadd_sim.sh 1 2;
# sh do_analysis.sh 1;

if [ "$#" -ne 1 ]; then
    echo "Error: Missing required parameter."
    echo "Usage: $0 <do_reweighting (0 or 1)>"
    exit 1
fi

do_reweighting=$1

# Check if the parameter is strictly 0 or 1
if [ "$do_reweighting" != "0" ] && [ "$do_reweighting" != "1" ]; then
    echo "Error: do_reweighting parameter must be exactly 0 or 1."
    echo "Usage: $0 <do_reweighting (0 or 1)>"
    exit 1
fi

# 1) Combine outputs
#echo "==> Running get_combinedoutput_mix.C"
root -l -q -b 'get_mixoutput.C(4,4)'
root -l -q -b 'get_combinedoutput_mix.C(4)'

# 2) Compute purity and efficiency
echo "==> Running get_purityefficiency.C"
#root -l -q -b 'get_purityefficiency.C(2)'
#root -l -q -b 'get_purityefficiency.C(3)'
root -l -q -b 'get_purityefficiency.C(4)'
#root -l -q -b 'get_purityefficiency.C(5)'
#root -l -q -b 'get_purityefficiency.C(6)'
#root -l -q -b 'get_purityefficiency.C(8)'

# 3) Perform unfolding iterations
echo "==> Running do_unfolding_iter.C"
#root -l -q -b 'do_unfolding_iter.C(2)'
#root -l -q -b 'do_unfolding_iter.C(3)'
root -l -q -b 'do_unfolding_iter.C(4)'
#root -l -q -b 'do_unfolding_iter.C(5)'
#root -l -q -b 'do_unfolding_iter.C(6)'
#root -l -q -b 'do_unfolding_iter.C(8)'

# (Optional) Get reweighting function for next iteration
if [ "$do_reweighting" -eq 0 ]; then
    echo "==> Running get_reweighthist.C (do_reweighting is 0)"
    #root -l -q -b 'get_reweighthist.C(2)'
    #root -l -q -b 'get_reweighthist.C(3)'
    root -l -q -b 'get_reweighthist.C(4)'
    #root -l -q -b 'get_reweighthist.C(5)'
    #root -l -q -b 'get_reweighthist.C(6)'
    #root -l -q -b 'get_reweighthist.C(8)'
else
    echo "==> Skipping get_reweighthist.C (do_reweighting is 1)"
fi

# 4) Generate the final spectrum
echo "==> Running get_finalspectrum.C"
root -l -q -b get_finalspectrum.C

# Done
echo "All analysis steps"

