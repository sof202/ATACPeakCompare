#!/bin/bash
#SBATCH --export=ALL
#SBATCH -p mrcq 
#SBATCH --time=12:00:00
#SBATCH -A Research_Project-MRC190311 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=16
#SBATCH --mem=10G
#SBATCH --mail-type=END 
#SBATCH --output=peakcompare_%j.log
#SBATCH --error=peakcompare_%j.err
#SBATCH --job-name=peakcompare

usage() {
cat << EOF
================================================================================
$(basename "$0")
================================================================================
Purpose: Compare peaks between two datasets over a specified region. Outputs
  a yet to be named metric. 
Arguments: \$1 -> full/relative path to config file
Author: Sam Fletcher
Contact: s.o.fletcher@exeter.ac.uk
================================================================================
EOF
    exit 0
}

validate_config_file() {
  local config_file_location
  config_file_location=$1
  local script_location
  script_location=$(\
    scontrol show job "${SLURM_JOB_ID}" | \
    grep "Command" | \
    awk '{print $1}' | \
    awk 'BEGIN {FS="="} {print $2}' \
  )
    PYTHON_SCRIPTS="$(dirname "${script_location}")/../Python_Scripts"
    PYTHON_SCRIPTS="$(realpath "${PYTHON_SCRIPTS}")"
  python3 \
    "${PYTHON_SCRIPTS}/validate_peak_compare_config.py" \
    "${config_file_location}"
  if [[ $? -eq 1 ]]; then
    echo "ERROR: Malformed config file detected."
    exit 1
  fi
}

move_log_files() {
  local log_directory
  log_directory="${LOG_DIRECTORY}/${USER}"
  local timestamp
  timestamp=$(date +%d-%h~%H-%M)
  mkdir -p "${log_directory}"
  mv "${SLURM_SUBMIT_DIR}/peakcompare_${SLURM_JOB_ID}.log" \
    "${log_directory}/${timestamp}_${SLURM_JOB_ID}_peakcompare.log"
  mv "${SLURM_SUBMIT_DIR}/peakcompare_${SLURM_JOB_ID}.err" \
    "${log_directory}/${timestamp}_${SLURM_JOB_ID}_peakcompare.err"
}

subset_file() {
  file=$1
  grep "${CHROMOSOME}" "${file}" > "${TEMP_DIRECTORY}/${file}"
}

main() {
    config_file=$1
    if [[ -f "${CONDA_EXE%/bin/conda}/etc/profile.d/conda.sh" ]]; then
        source "${CONDA_EXE%/bin/conda}/etc/profile.d/conda.sh"
    else
        # fallback
        source "${CONDA_EXE%/condabin/conda}/etc/profile.d/conda.sh" || \
            { >&2 echo "Could not find conda executable, please check you have activated conda first."
            exit 1; }
    fi
    conda activate PeakCompare-MACS || \
        { >&2 echo "Conda environment does not exist. Please run setup script first."
        exit 1; }
    validate_config_file "${config_file}"
    source "${config_file}" || exit 1
    move_log_files
    mkdir -p "${TEMP_DIRECTORY}"

    subset_file "${REFERENCE_MERGED_PEAK_FILE}"
    subset_file "${REFERENCE_UNMERGED_PEAK_FILE}"
    subset_file "${REFERENCE_BIAS_TRACK_FILE}"
    subset_file "${REFERENCE_COVERAGE_TRACK_FILE}"
    subset_file "${COMPARISON_BIAS_TRACK_FILE}"
    subset_file "${COMPARISON_COVERAGE_TRACK_FILE}"
    subset_file "${COMPARISON_PVALUE_FILE}"

    python3 \
        "${PYTHON_SCRIPTS}/peak_compare.py" \
        "$(if [[ ${UNMERGED} -eq 1 ]]; then echo "--unmerged"; fi)" \
        "${CHROMOSOME}" \
        "${START}" \
        "${END}" \
        "${TEMP_DIRECTORY}/${REFERENCE_MERGED_PEAK_FILE}" \
        "${TEMP_DIRECTORY}/${REFERENCE_UNMERGED_PEAK_FILE}" \
        "${TEMP_DIRECTORY}/${REFERENCE_BIAS_TRACK_FILE}" \
        "${TEMP_DIRECTORY}/${REFERENCE_COVERAGE_TRACK_FILE}" \
        "${TEMP_DIRECTORY}/${COMPARISON_BIAS_TRACK_FILE}" \
        "${TEMP_DIRECTORY}/${COMPARISON_COVERAGE_TRACK_FILE}" \
        "${TEMP_DIRECTORY}/${COMPARISON_PVALUE_FILE}" \
        "${CUTOFF}"

    rm "${TEMP_DIRECTORY}/${REFERENCE_SAMPLE_NAME}*.bed" \
      "${TEMP_DIRECTORY}/${COMPARISON_SAMPLE_NAME}*.bed"
}

if [[ $# -ne 1 ]]; then usage; fi
main "$1"
