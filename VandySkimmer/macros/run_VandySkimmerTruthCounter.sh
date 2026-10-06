#! /bin/bash

caloDSTlist=${1:-""}
truthDSTlist=${2:-""}
jetDSTlist=${3:-""}
globalDSTlist=${4:-""}
outDir=${5:-""}
installdir=${6:-$MYINSTALL}
thisDir=${7:-""}
Pythia_or_Herwig=${8:-'true'}
source /opt/sphenix/core/bin/sphenix_setup.sh -n ana.542
source /opt/sphenix/core/bin/setup_local.sh $installdir

n_files=`wc -l < $caloDSTlist`
root -x -q ${thisDir}/Fun4All_VandySkimmerTruthCounter.C\(\"$caloDSTlist\",\"$jetDSTlist\",\"$truthDSTlist\",\"$globalDSTlist\",\"$outDir\",\"$n_files\",\"$Pythia_or_Herwig\"\)

