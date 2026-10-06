#!/usr/bin/env bash

# Exit on errors, undefined variables, and pipeline failures so a bad
# submission state does not silently continue to the next scan tag.
set -euo pipefail

analysis_dir="/sphenix/user/hjheng/sPHENIXRepo/analysis/ChargedHadronMultiplicity-SmallSystem/analysis"
job_file="${analysis_dir}/condor_TrackletAnaCombinatoric.job"
scan_root="/sphenix/tg/tg01/hf/hjheng/ppg-dNdEta-OOpp/simulation-nopileup-ntuple/ACTS"
poll_interval=60
dry_run=0

if [[ "${1:-}" == "--dry-run" ]]; then
    dry_run=1
fi

if [[ ! -f "${job_file}" ]]; then
    echo "Missing job file: ${job_file}" >&2
    exit 1
fi

if [[ ! -d "${scan_root}" ]]; then
    echo "Missing scan root: ${scan_root}" >&2
    exit 1
fi

original_scan_tag=$(
    awk '
        $1 == "scanTag" && $2 == "=" {
            print $3
            exit
        }
    ' "${job_file}"
)

if [[ -z "${original_scan_tag}" ]]; then
    echo "Could not read the current scanTag from ${job_file}" >&2
    exit 1
fi

restore_original_scantag() {
    sed -i "s|^scanTag[[:space:]]*=.*$|scanTag            = ${original_scan_tag}|" "${job_file}"
}

# Always restore the original scanTag when the script exits, even if it
# stops early because submission or queue monitoring failed.
trap restore_original_scantag EXIT

# Read all immediate subdirectory names under the ACTS scan area into a Bash
# array. `-printf '%f\n'` keeps only the basename, and `sort` makes the run
# order deterministic.
mapfile -t scan_tags < <(find "${scan_root}" -mindepth 1 -maxdepth 1 -type d -printf '%f\n' | sort)

if [[ ${#scan_tags[@]} -eq 0 ]]; then
    echo "No scan tags found under ${scan_root}" >&2
    exit 1
fi

wait_for_cluster() {
    local cluster_id="$1"

    while true; do
        local running_jobs held_jobs

        running_jobs=$(
            condor_q -constraint "ClusterId == ${cluster_id}" -autoformat ClusterId ProcId 2>/dev/null \
                | awk 'NF {count++} END {print count + 0}'
        )

        held_jobs=$(
            condor_q -hold -constraint "ClusterId == ${cluster_id}" -autoformat ClusterId ProcId 2>/dev/null \
                | awk 'NF {count++} END {print count + 0}'
        )

        if [[ "${held_jobs}" -gt 0 ]]; then
            echo "Cluster ${cluster_id} has ${held_jobs} held job(s)." >&2
            echo "Stopping scan-tag submission loop so you can inspect the held jobs." >&2
            exit 1
        fi

        if [[ "${running_jobs}" -eq 0 ]]; then
            break
        fi

        echo "Cluster ${cluster_id}: ${running_jobs} job(s) still in queue. Sleeping ${poll_interval}s."
        sleep "${poll_interval}"
    done
}

cd "${analysis_dir}"

for scan_tag in "${scan_tags[@]}"; do
    echo "Submitting scanTag=${scan_tag}"

    sed -i "s|^scanTag[[:space:]]*=.*$|scanTag            = ${scan_tag}|" "${job_file}"

    if [[ "${dry_run}" -eq 1 ]]; then
        echo "[DRY-RUN] Would run: condor_submit ${job_file}"
        echo "[DRY-RUN] Current scanTag line:"
        awk '$1 == "scanTag" && $2 == "=" { print }' "${job_file}"
        continue
    fi

    submit_output=$(condor_submit "${job_file}")
    echo "${submit_output}"

    cluster_id=$(
        awk '
            match($0, /cluster[[:space:]]+([0-9]+)/, m) {
                print m[1]
                exit
            }
        ' <<< "${submit_output}"
    )

    if [[ -z "${cluster_id}" ]]; then
        echo "Failed to parse cluster id from condor_submit output." >&2
        exit 1
    fi

    wait_for_cluster "${cluster_id}"
    echo "Cluster ${cluster_id} finished for scanTag=${scan_tag}"
done

echo "All scan tags completed."
