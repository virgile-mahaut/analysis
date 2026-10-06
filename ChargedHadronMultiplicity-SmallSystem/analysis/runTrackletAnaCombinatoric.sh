#! /bin/bash
export USER="$(id -u -n)"
export LOGNAME=${USER}
export HOME=/sphenix/u/${LOGNAME}
source /opt/sphenix/core/bin/sphenix_setup.sh -n new

export MYINSTALL=$HOME/install
export LD_LIBRARY_PATH=$MYINSTALL/lib:$LD_LIBRARY_PATH
# export ROOT_INCLUDE_PATH=$MYINSTALL/include:$ROOT_INCLUDE_PATH

source /opt/sphenix/core/bin/setup_local.sh $MYINSTALL

# This is has to be after sourcing the setup_local.sh, otherwise ROOT_INCLUDE_PATH will be overwritten by the setup scripts
# export WORKSPACE=/sphenix/user/hjheng/sPHENIXRepo
# export ROOT_INCLUDE_PATH=$WORKSPACE/macros/common:$ROOT_INCLUDE_PATH

# print the environment - needed for debugging
# printenv
# $3 is the output filename, get the the directory of the output file, and make the directory if it doesn't exist
outdir=$(dirname "$3")
if [ ! -d "$outdir" ]; then
    mkdir -p "$outdir"
fi
./TrackletAna_Combinatoric $1 $2 $3 $4 $5 $6 $7 $8

echo all done
