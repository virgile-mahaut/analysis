#!/usr/bin/env bash

# Exit on errors, undefined variables, and pipeline failures so a failed
# merge or ROOT macro call stops the scan loop immediately.
set -euo pipefail

analysis_dir="/sphenix/user/hjheng/sPHENIXRepo/analysis/ChargedHadronMultiplicity-SmallSystem/analysis"
hist_root="/sphenix/tg/tg01/hf/hjheng/ppg-dNdEta-OOpp/TrackletAna-HistOutput/simulation/ACTS"
macro_file="${analysis_dir}/plot_tklCombinatoric.C"
dry_run=0

if [[ "${1:-}" == "--dry-run" ]]; then
    dry_run=1
fi

if [[ ! -d "${hist_root}" ]]; then
    echo "Missing histogram root directory: ${hist_root}" >&2
    exit 1
fi

if [[ ! -f "${macro_file}" ]]; then
    echo "Missing macro file: ${macro_file}" >&2
    exit 1
fi

# Read all immediate scan-tag subdirectories and process them in a stable order.
mapfile -t scan_tags < <(find "${hist_root}" -mindepth 1 -maxdepth 1 -type d -printf '%f\n' | sort)

if [[ ${#scan_tags[@]} -eq 0 ]]; then
    echo "No scan-tag directories found under ${hist_root}" >&2
    exit 1
fi

needs_merge() {
    local scan_dir="$1"
    local merged_file="${scan_dir}/merged_hist.root"

    shopt -s nullglob
    local hist_files=("${scan_dir}"/hist_*.root)
    shopt -u nullglob

    if [[ ${#hist_files[@]} -eq 0 ]]; then
        echo "No individual histogram files found in ${scan_dir}" >&2
        return 2
    fi

    if [[ ! -f "${merged_file}" ]]; then
        return 0
    fi

    # Re-merge if any individual histogram file is newer than merged_hist.root.
    if find "${scan_dir}" -maxdepth 1 -name 'hist_*.root' -newer "${merged_file}" -print -quit | grep -q .; then
        return 0
    fi

    return 1
}

run_plot_macro() {
    local scan_tag="$1"
    local merged_file="$2"
    local outstem="plot_tklCombinatoric_${scan_tag}"
    local plotdir="${analysis_dir}/figure/figure-tklCombinatoric-simulation-${scan_tag}"
    local root_cmd

    root_cmd=$(printf '%s("%s","tkl_Combinatoric","absdPhi",true,-10.0,10.0,0.0,0.015,0.0,0.2,"%s","%s",true,"colz text",".3f",true)' \
        "${macro_file}" "${merged_file}" "${outstem}" "${plotdir}")

    if [[ "${dry_run}" -eq 1 ]]; then
        echo "[DRY-RUN] Would run ROOT macro for scanTag=${scan_tag}"
        echo "[DRY-RUN] root -l -b -q '${root_cmd}'"
        return 0
    fi

    echo "Running plot_tklCombinatoric for scanTag=${scan_tag}"
    root -l -b -q "${root_cmd}"
}

cd "${analysis_dir}"

for scan_tag in "${scan_tags[@]}"; do
    scan_dir="${hist_root}/${scan_tag}"
    merged_file="${scan_dir}/merged_hist.root"

    echo "Processing scanTag=${scan_tag}"

    merge_status=0
    if needs_merge "${scan_dir}"; then
        merge_status=0
    else
        merge_status=$?
    fi

    if [[ ${merge_status} -eq 2 ]]; then
        echo "Skipping scanTag=${scan_tag} because no hist_*.root files were found." >&2
        continue
    fi

    if [[ ${merge_status} -eq 0 ]]; then
        if [[ "${dry_run}" -eq 1 ]]; then
            echo "[DRY-RUN] Would run: (cd ${scan_dir} && hadd-toplevel -f -k -j 20 merged_hist.root hist_*.root)"
        else
            echo "Merging histograms in ${scan_dir}"
            (
                cd "${scan_dir}"
                hadd-toplevel -f -k -j 20 merged_hist.root hist_*.root
            )
        fi
    else
        echo "merged_hist.root is already up to date for scanTag=${scan_tag}"
    fi

    run_plot_macro "${scan_tag}" "${merged_file}"
done

echo "Finished processing all scan tags."
