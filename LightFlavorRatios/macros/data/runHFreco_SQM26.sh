#!/bin/bash

#source /opt/sphenix/core/bin/sphenix_setup.sh -n new 

#export USER="$(id -u -n)"
#export LOGNAME=${USER}
#export HOME=/sphenix/u/${LOGNAME}

#export SPHENIX=${HOME}/sPHENIX
#export MYINSTALL=$SPHENIX/install
#export LD_LIBRARY_PATH=$MYINSTALL/lib:$LD_LIBRARY_PATH
#export ROOT_INCLUDE_PATH=$MYINSTALL/include:$ROOT_INCLUDE_PATH
#
#source /opt/sphenix/core/bin/setup_local.sh $MYINSTALL

nEvents=$1
inDst=$2
outDir=$3
nSkip=$4

#if [[ "${inDst}" == *.root ]]; then
#  getinputfiles.pl $inDst
#elif [[ "${inDst}" == *.list ]]; then
#  getinputfiles.pl --filelist $inDst
#fi

# print the environment - needed for debugging
#printenv

echo running: runHFreco_SQM26.sh $*
root.exe -q -b Fun4All_HF_SQM26.C\(${nEvents},\"${inDst}\",\"${outDir}\",${nSkip}\)
echo Script done