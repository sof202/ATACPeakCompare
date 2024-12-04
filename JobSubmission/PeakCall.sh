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

move_log_files() {
  log_directory="${OUTPUT_DIRECTORY}/LogFiles/${USER}"
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
}

build_model() {
    macs3 predictd \
        -i "${OUTPUT_DIRECTORY}/${SAMPLE_NAME}_filtered.bed" \
        -f "${FILE_TYPE}" \
        -g "${GENOME_SIZE}" \
        -m "${MFOLD_LOWER}" "${MFOLD_UPPER}" \
        --outdir "${OUTPUT_DIRECTORY}" 2> \
        "${OUTPUT_DIRECTORY}/model.txt"

    if grep -qi "mfold" "${OUTPUT_DIRECTORY}/model.txt"; then
cat 1>&2 << EOF
WARNING
MACS was unable to build the model.
Consider increasing the range of MFOLD variables in config file.
EOF
    exit 1
    fi
}

get_fragment_length() {
    fragment_length=$(\
        grep "predicted fragment" "${OUTPUT_DIRECTORY}/model.txt" | \
        grep -Po "\d+ bps" | \
        grep -Po "\d+" \
    )
}

get_read_length() {
    read_length=$(\
        grep "tag size" "${OUTPUT_DIRECTORY}/model.txt" | \
        grep -Po "\d+ bps" | \
        grep -Po "\d+" \
    )
}

get_number_of_reads() {
    number_of_reads=$(\
        grep "total reads" "${OUTPUT_DIRECTORY}/model.txt" | \
        grep -Po ": \d+" | \
        grep -Po "\d+" \
    )
}

get_coverage_track() {
    macs3 pileup \
        -i "${OUTPUT_DIRECTORY}/${SAMPLE_NAME}_filtered.bed" \
        -f "${FILE_TYPE}" \
        --extsize "${fragment_length}" \
        -o "${OUTPUT_DIRECTORY}/${SAMPLE_NAME}_pileup.bdg"
}

get_bias_track() {
    background_file=$1
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
        -o "${OUTPUT_DIRECTORY}/bias_track.bdg"

    if [[ "${DEBUG_MODE}" -ne 1 ]]; then
        rm "${OUTPUT_DIRECTORY}"/background*.bdg
    fi
}

main() {
    config_file=$1
    validate_config_file "${config_file}"
    source "${config_file}" || exit 1
    move_log_files
    if [[ -f "${CONDA_EXE%/bin/conda}/etc/profile.d/conda.sh" ]]; then
        source "${CONDA_EXE%/bin/conda}/etc/profile.d/conda.sh"
    else
        # fallback
        source "${CONDA_EXE%/condabin/conda}/etc/profile.d/conda.sh" || \
            { >&2 echo "Could not find conda executable, please check you have activated conda first."
            exit 1; }
    fi
    conda activate MACS3
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
    fi
    get_coverage_track
    if [[ -f "${CONTROL_FILE}" ]]; then
        get_bias_track "${CONTROL_FILE}"
        scale_bias_track
    else
        get_bias_track "${INPUT_FILE}"
    fi
    get_p_values
    get_cutoff_analysis
    get_best_cutoff
    get_unmerged_peaks
    get_merged_peaks
}

if [[ $#1 -ne 1 ]]; then usage; fi
main "$1"
