#!/bin/bash
#
#

if [[ $# -lt 1 ]]
then
  echo 'Usage: runupcmeson.cmd DST_TRKR_TRACKS.root <nevents>'
  exit -1
fi

ulimit -c 0	# no core files

echo PWD=${PWD}
echo LD_LIBRARY_PATH=${LD_LIBRARY_PATH}
echo HOST=`hostname`

TOPDIR=${PWD}

# get dst or list fname
fname=$1

nevents=0
if [[ $# -eq 2 ]]
then
  nevents=$2
fi

# output save directory
if [[ ${fname} =~ .root$ ]]
then
  runno=$(getrunfromfname.sh ${fname})
  seedf=${fname}
  clusf=$( echo ${fname} | sed 's/SEED/CLUSTER/' )
elif [[ ${fname} =~ .list$ ]]
then
  runno=$(getrunfromfname.sh $(head -1 ${fname}))
fi

SAVEDIR=${TOPDIR}/OUTPUT/${runno}/${fname%.root}
mkdir -p ${SAVEDIR}
cd ${SAVEDIR}
ln -sf ${TOPDIR}/Fun4All_UPCMeson.C .
ln -sf ${TOPDIR}/Fun4All_TrackFitting.C .
ln -sf ${TOPDIR}/my_Fun4All_TrackFitting.C .
ln -sf ${TOPDIR}/runupcmeson.cmd .

if [[ ! -z ${_CONDOR_SCRATCH_DIR} ]]
then

  mkdir -p ${_CONDOR_SCRATCH_DIR}
  cp -p *.C ${_CONDOR_SCRATCH_DIR}
  cd ${_CONDOR_SCRATCH_DIR}

  if [[ ${fname} =~ .root$ ]]
  then
    getinputfiles.pl ${seedf}
    getinputfiles.pl ${clusf}
  elif [[ ${fname} =~ .list$ ]]
  then
    getinputfiles.pl --filelist ${fname}
    sed 's/SEED/CLUSTER/' ${fname} > /tmp/clus.$$
    getinputfiles.pl --filelist /tmp/clus.$$
  else
    echo "Unknown file type ${fname}, exiting"
    exit 1
  fi

fi

if [[ ${fname##*/} =~ ^dst_tracks || ${fname##*/} =~ ^DST_TRKR_TRACKS ]]
then
  ln -sf ${TOPDIR}/${fname} .
  echo root.exe -b -q Fun4All_UPCMeson.C\(${nevents},\"${fname}\"\)
  root.exe -b -q Fun4All_UPCMeson.C\(${nevents},\"${fname}\"\)
elif [[ ${fname} =~ DST_TRKR_SEED ]]
then
  outfname="aaa"
  echo root.exe -b -q my_Fun4All_TrackFitting.C\($nevents,\"${seedf}\",\"${clusf}\",\"${outfname}\",true\)
  root.exe -b -q my_Fun4All_TrackFitting.C\(${nevents},\"${seedf}\",\"${clusf}\",\"${outfname}\",true\)
fi


if [[ ! -z "${_CONDOR_SCRATCH_DIR}" ]]
then
  cd ${SAVEDIR}
  cp -p ${_CONDOR_SCRATCH_DIR}/upcmeson*.root .
  cp -p ${_CONDOR_SCRATCH_DIR}/aaa*.root .
fi

