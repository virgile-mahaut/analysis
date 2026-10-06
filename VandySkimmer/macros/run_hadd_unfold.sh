#!/bin/bash                                                                                                                                                                                               

export USER="$(id -u -n)"
export LOGNAME=${USER}
export HOME=/sphenix/u/${LOGNAME}/macros/detectors/sPHENIX/
export MYINSTALL=/sphenix/user/bkimelman/sPHENIX/install/

source /opt/sphenix/core/bin/sphenix_setup.sh -n
source /opt/sphenix/core/bin/setup_local.sh $MYINSTALL


set -euo pipefail
shopt -s nullglob

baseDir="$2"
fileBase="$3"
closureMode="$4"

closureModeStr=""
if [[ "$closureMode" == "kFull" ]]; then
    closureModeStr="fullClosure"
elif [[ "$closureMode" == "kHalf" ]]; then
    closureModeStr="halfClosure"
elif [[ "$closureMode" == "kNoTest" ]]; then
    closureModeStr="noTest"
fi

jetSample=$1
chunkSize=10

rawDir="$baseDir/Jet${jetSample}"
workDir="$baseDir/.merge_tmp_Jet${jetSample}-${closureModeStr}"
finalOut="$baseDir/${fileBase}_Jet${jetSample}-${closureModeStr}.root"
#finalOut="$baseDir/${fileBase}_Jet${jetSample}-all.root"

mkdir -p "$workDir"
rm -f -- "$workDir"/*.root

extract_range() {
    local base="$1"

    # Raw files
    if [[ $base =~ ^${fileBase}_Jet${jetSample}_seg([0-9]{6})_to_([0-9]{6})-${closureModeStr}\.root$ ]];
    #if [[ $base =~ ^${fileBase}_Jet${jetSample}_seg([0-9]{6})_to_([0-9]{6})\.root$ ]];
    then
	printf '%s %s\n' "${BASH_REMATCH[1]}" "${BASH_REMATCH[2]}"
	return 0
    fi

    # Intermediate files
    if [[ $base =~ ^merged_Jet${jetSample}_([0-9]{6})_to_([0-9]{6})-${closureModeStr}\.root$ ]];
    #if [[ $base =~ ^merged_Jet${jetSample}_([0-9]{6})_to_([0-9]{6})\.root$ ]];
    then
	printf '%s %s\n' "${BASH_REMATCH[1]}" "${BASH_REMATCH[2]}"
	return 0
    fi

    echo "ERROR: cannot parse range from: $base" >&2
    return 1
}

current=( "$rawDir/${fileBase}_Jet${jetSample}_"*"-${closureModeStr}.root" )
#current=( "$rawDir/${fileBase}_Jet${jetSample}_"*".root" )

if (( ${#current[@]} == 0 ));
then
    echo "ERROR: no input files found" >&2
    exit 1
fi

delete_current=false
round=0

while (( ${#current[@]} > 1)); do
    echo "Round $round: ${#current[@]} inputs"

    next=()

    for ((i=0; i<${#current[@]}; i+=chunkSize)); do
	chunk=( "${current[@]:i:chunkSize}" )

	firstBase=$(basename "${chunk[0]}")
        lastIndex=$((${#chunk[@]} - 1))
        lastBase=$(basename "${chunk[$lastIndex]}")

        read -r firstNum _ < <(extract_range "$firstBase")
        read -r _ secondNum < <(extract_range "$lastBase")
	
	outfile="$workDir/merged_Jet${jetSample}_${firstNum}_to_${secondNum}-${closureModeStr}.root"
	#outfile="$workDir/merged_Jet${jetSample}_${firstNum}_to_${secondNum}.root"

	echo "   Round $round: hadd -f $outfile"
	hadd -f "$outfile" "${chunk[@]}"

	next+=( "$outfile" )
    done

    echo "Done with round $round"
    
    printf 'Next round inputs:\n%s\n' "${next[@]}"

    if (( round > 0 )); then
        rm -f -- "${current[@]}"
    fi

    echo "Round $round produced ${#next[@]} merged files"
    current=( "${next[@]}" )
    echo "After round $round, current has ${#current[@]} files"
    ((++round))
done

mv -f -- "${current[0]}" "$finalOut"
rm -rf -- "$workDir"

echo "Final file: $finalOut"
