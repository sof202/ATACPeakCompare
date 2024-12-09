#!/bin/bash
#SBATCH --export=ALL
#SBATCH -p mrcq 
#SBATCH --time=12:00:00
#SBATCH -A Research_Project-MRC190311 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=16
#SBATCH --mem=10G
#SBATCH --mail-type=END 
#SBATCH --output=peakcall_%j.log
#SBATCH --error=peakcall_%j.err
#SBATCH --job-name=peakcall

usage() {
cat << EOF
================================================================================
$(basename "$0")
================================================================================
Purpose: Build MACS3 model and then obtain 4 bedgraph files:
  - The bias track (lambdas)
  - p-values for each basepair
  - Unmerged peaks
  - Merged peaks
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
    "${PYTHON_SCRIPTS}/validate_config_file.py" \
    "${config_file_location}"
  if [[ $? -eq 1 ]]; then
    echo "ERROR: Malformed config file detected."
    exit 1
  fi
}

move_log_files() {
  local log_directory
  log_directory="${OUTPUT_DIRECTORY}/LogFiles/${USER}"
  local timestamp
  timestamp=$(date +%d-%h~%H-%M)
  mkdir -p "${log_directory}"
  mv "${SLURM_SUBMIT_DIR}/peakcall_${SLURM_JOB_ID}.log" \
    "${log_directory}/${timestamp}_${SLURM_JOB_ID}_peakcall.log"
  mv "${SLURM_SUBMIT_DIR}/peakcall_${SLURM_JOB_ID}.err" \
    "${log_directory}/${timestamp}_${SLURM_JOB_ID}_peakcall.err"
}

remove_duplicates() {
    macs3 filterdup \
        -i "${INPUT_FILE}" \
        -f "${FILE_TYPE}" \
        -o "${OUTPUT_DIRECTORY}/${SAMPLE_NAME}_filtered.bed"
    if [[ -f "${CONTROL_FILE}" ]]; then
        macs3 filterdup \
            -i "${CONTROL_FILE}" \
            -f "${FILE_TYPE}" \
            -o "${OUTPUT_DIRECTORY}/${SAMPLE_NAME}_control_filtered.bed"
    fi
}

build_model() {
    macs3 predictd \
        -i "${OUTPUT_DIRECTORY}/${SAMPLE_NAME}_filtered.bed" \
        -f "${FILE_TYPE}" \
        -g "${GENOME_SIZE}" \
        -m "${MFOLD_LOWER}" "${MFOLD_UPPER}" \
        --outdir "${OUTPUT_DIRECTORY}" 2> \
        "${OUTPUT_DIRECTORY}/${SAMPLE_NAME}_model.txt"

    if [[ -f "${CONTROL_FILE}" ]]; then
        macs3 predictd \
            -i "${OUTPUT_DIRECTORY}/${SAMPLE_NAME}_control_filtered.bed" \
            -f "${FILE_TYPE}" \
            -g "${GENOME_SIZE}" \
            -m "${MFOLD_LOWER}" "${MFOLD_UPPER}" \
            --outdir "${OUTPUT_DIRECTORY}" 2> \
            "${OUTPUT_DIRECTORY}/${SAMPLE_NAME}_control_model.txt"
    fi

if grep -qi "mfold" "${OUTPUT_DIRECTORY}"/*model.txt; then
cat 1>&2 << EOF
WARNING
MACS was unable to build the model.
Consider increasing the range of MFOLD variables in config file.
Alternatively, you can specify your own model parameters in the config file.
If you do this, be sure to set \$BUILD_MODEL to 0.
EOF
exit 1
fi
}

# These functions are very messy, but there is not easily parsable model that
# MACS can create (unless you modify the source code it would appear).
get_fragment_length() {
    fragment_length=$(\
        grep "predicted fragment" "${OUTPUT_DIRECTORY}/${SAMPLE_NAME}_model.txt" | \
        grep -Po "\d+ bps" | \
        grep -Po "\d+" \
    )
    if [[ -f "${CONTROL_FILE}" ]]; then
        control_fragment_length=$(\
            grep "predicted fragment" "${OUTPUT_DIRECTORY}/${SAMPLE_NAME}_control_model.txt" | \
            grep -Po "\d+ bps" | \
            grep -Po "\d+" \
        )
    fi
}

get_read_length() {
    read_length=$(\
        grep "tag size" "${OUTPUT_DIRECTORY}/${SAMPLE_NAME}_model.txt" | \
        grep -Po "\d+ bps" | \
        grep -Po "\d+" \
    )
}

get_number_of_reads() {
    number_of_reads=$(\
        grep "total reads" "${OUTPUT_DIRECTORY}/${SAMPLE_NAME}_model.txt" | \
        grep -Po ": \d+" | \
        grep -Po "\d+" \
    )
    if [[ -f "${CONTROL_FILE}" ]]; then
        number_of_control_reads=$(\
            grep "total reads" "${OUTPUT_DIRECTORY}/${SAMPLE_NAME}_control_model.txt" | \
            grep -Po ": \d+" | \
            grep -Po "\d+" \
        )
    fi
}

get_coverage_track() {
    macs3 pileup \
        -i "${OUTPUT_DIRECTORY}/${SAMPLE_NAME}_filtered.bed" \
        -f "${FILE_TYPE}" \
        --extsize "${fragment_length}" \
        -o "${OUTPUT_DIRECTORY}/${SAMPLE_NAME}_coverage.bdg"
}

get_bias_track() {
    background_file=$1
    local fragment_length
    fragment_length=$2
    local number_of_reads
    number_of_reads=$3

    macs3 pileup \
        -i "${background_file}" \
        -f "${FILE_TYPE}" \
        -B \
        --extsize "$(echo "scale=5; ${fragment_length}/2" | bc)" \
        -o "${OUTPUT_DIRECTORY}/background_fragment.bdg"
    macs3 pileup \
        -i "${background_file}" \
        -f "${FILE_TYPE}" \
        -B \
        --extsize "$(echo "scale=5; ${SMALL_LOCAL_SIZE}/2" | bc)" \
        -o "${OUTPUT_DIRECTORY}/background_small_local.bdg"
    macs3 pileup \
        -i "${background_file}" \
        -f "${FILE_TYPE}" \
        -B \
        --extsize "$(echo "scale=5; ${LARGE_LOCAL_SIZE}/2" | bc)" \
        -o "${OUTPUT_DIRECTORY}/background_large_local.bdg"

    macs3 bdgopt \
        -i "${OUTPUT_DIRECTORY}/background_small_local.bdg" \
        -m multiply \
        -p "$(echo "scale=5; ${fragment_length}/${SMALL_LOCAL_SIZE}" | bc)" \
        -o "${OUTPUT_DIRECTORY}/background_small_local_norm.bdg"
    macs3 bdgopt \
        -i "${OUTPUT_DIRECTORY}/background_large_local.bdg" \
        -m multiply \
        -p "$(echo "scale=5; ${fragment_length}/${LARGE_LOCAL_SIZE}" | bc)" \
        -o "${OUTPUT_DIRECTORY}/background_large_local_norm.bdg"

    macs3 bdgcmp \
        -m max \
        -t "${OUTPUT_DIRECTORY}/background_small_local_norm.bdg" \
        -c "${OUTPUT_DIRECTORY}/background_large_local_norm.bdg" \
        -o "${OUTPUT_DIRECTORY}/background_small_large_local.bdg"
    macs3 bdgcmp \
        -m max \
        -t "${OUTPUT_DIRECTORY}/background_small_large_local.bdg" \
        -c "${OUTPUT_DIRECTORY}/background_fragment.bdg" \
        -o "${OUTPUT_DIRECTORY}/background_local.bdg"

    global_background="$(\
        echo "scale=5; ${number_of_reads}*${fragment_length}/${GENOME_SIZE}" | \
        bc \
    )"
    macs3 bdgopt \
        -i "${OUTPUT_DIRECTORY}/background_local.bdg" \
        -m max \
        -p "${global_background}" \
        -o "${OUTPUT_DIRECTORY}/${SAMPLE_NAME}_bias_track.bdg"

    if [[ "${DEBUG_MODE}" -ne 1 ]]; then
        rm "${OUTPUT_DIRECTORY}"/background*.bdg
    fi
}

scale_bias_track() {
    read_ratio=$(\
        echo "scale=5; ${number_of_reads}/${number_of_control_reads}" | \
        bc \
    )
    macs3 bdgopt \
        -i "${OUTPUT_DIRECTORY}/${SAMPLE_NAME}_bias_track.bdg" \
        -m multiply \
        -p "${read_ratio}" \
        -o "${OUTPUT_DIRECTORY}/${SAMPLE_NAME}_bias_track_scaled.bdg"
    rm "${OUTPUT_DIRECTORY}/${SAMPLE_NAME}_bias_track.bdg"
    mv \
        "${OUTPUT_DIRECTORY}/${SAMPLE_NAME}_bias_track_scaled.bdg" \
        "${OUTPUT_DIRECTORY}/${SAMPLE_NAME}_bias_track.bdg" 
}

get_p_values() {
    macs3 bdgcmp \
        -t "${OUTPUT_DIRECTORY}/${SAMPLE_NAME}_coverage.bdg" \
        -c "${OUTPUT_DIRECTORY}/bias_track.bdg" \
        -m ppois \
        -o "${OUTPUT_DIRECTORY}/${SAMPLE_NAME}_pvalues.bdg"
}

get_cutoff_analysis() {
    macs3 bdgpeakcall \
        -i "${OUTPUT_DIRECTORY}/${SAMPLE_NAME}_pvalues.bdg" \
        --cutoff-analysis \
        -o "${OUTPUT_DIRECTORY}/${SAMPLE_NAME}_cutoff_analysis.txt"
}

get_best_cutoff() {
    # We need to use tail here in case of any warning messages due to dtypes
    # not being configured correctly for pandas
    cutoff=$(\
        python3 "${PYTHON_SCRIPTS}/get_cutoff.py" \
        "${OUTPUT_DIRECTORY}/${SAMPLE_NAME}_cutoff_analysis.txt" \
        "${AVERAGE_PEAK_LENGTH}" | \
        tail -1 \
    )
    echo "Using a cutoff value of ${cutoff}."
}

get_unmerged_peaks() {
    macs3 bdgpeakcall \
        -i "${OUTPUT_DIRECTORY}/${SAMPLE_NAME}_pvalues.bdg" \
        -c "${cutoff}" \
        -l "${fragment_length}" \
        -o "${OUTPUT_DIRECTORY}/${SAMPLE_NAME}_unmerged_peaks.bed"
}

get_merged_peaks() {
    macs3 bdgpeakcall \
        -i "${OUTPUT_DIRECTORY}/${SAMPLE_NAME}_pvalues.bdg" \
        -c "${cutoff}" \
        -l "${fragment_length}" \
        -g "${read_length}" \
        -o "${OUTPUT_DIRECTORY}/${SAMPLE_NAME}_merged_peaks.bed"
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
    remove_duplicates
    if [[ "${BUILD_MODEL}" -eq 1 ]]; then
        build_model
        get_fragment_length 
        get_read_length
        get_number_of_reads
    else
        fragment_length="${FRAGMENT_LENGTH}"
        read_length="${READ_LENGTH}"
        number_of_reads="${NUMBER_OF_READS}"
        control_fragment_length="${CONTROL_FRAGMENT_LENGTH}"
        number_of_control_reads="${NUMBER_OF_CONTROL_READS}"
    fi
    get_coverage_track
    if [[ -f "${CONTROL_FILE}" ]]; then
        get_bias_track \
            "${CONTROL_FILE}" \
            "${control_fragment_length}" \
            "${number_of_control_reads}"
        scale_bias_track
    else
        get_bias_track \
            "${INPUT_FILE}" \
            "${fragment_length}" \
            "${number_of_reads}"
    fi
    get_p_values
    if [[ -z "${CUTOFF}" ]]; then
        get_cutoff_analysis
        get_best_cutoff
    else
        cutoff="${CUTOFF}"
    fi
    get_unmerged_peaks
    get_merged_peaks
}

if [[ $#1 -ne 1 ]]; then usage; fi
main "$1"
