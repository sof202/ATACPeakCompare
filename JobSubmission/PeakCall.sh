#!/bin/bash
#SBATCH --export=ALL
#SBATCH -p mrcq 
#SBATCH --time=12:00:00
#SBATCH -A Research_Project-MRC190311 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=16
#SBATCH --mem=10G
#SBATCH --mail-type=END 
#SBATCH --output=peakcall%j.log
#SBATCH --error=peakcall%j.err
#SBATCH --job-name=peakcall

remove_duplicates() {
    macs3 filterdup \
        -i "${INPUT_SAMPLE}" \
        -f "${FILE_TYPE}" \
        -o "${PROCESSING_DIRECTORY}/${SAMPLE_NAME}_filtered.bed"
}

build_model() {
    macs3 predictd \
        -i "${PROCESSING_DIRECTORY}/${SAMPLE_NAME}_filtered.bed" \
        -f "${FILE_TYPE}" \
        -g "${GENOME_SIZE}" \
        -m "${MFOLD_LOWER}" "${MFOLD_UPPER}" \
        --outdir "${PROCESSING_DIRECTORY}" 2> \
        "${PROCESSING_DIRECTORY}/model.txt"

    if grep -qi "mfold" "${PROCESSING_DIRECTORY}/model.txt"; then
cat >> "${LOG_FILE}" << EOF
WARNING
MACS was unable to build the model.
Consider increasing the range of MFOLD variables in config file.
EOF
    exit 1
    fi
}

get_fragment_length() {
    fragment_length=$(\
        grep "predicted fragment" "${PROCESSING_DIRECTORY}/model.txt" | \
        grep -Po "\d+ bps" | \
        grep -Po "\d+" \
    )
}

get_read_length() {
    read_length=$(\
        grep "tag size" "${PROCESSING_DIRECTORY}/model.txt" | \
        grep -Po "\d+ bps" | \
        grep -Po "\d+" \
    )
}

get_number_of_reads() {
    number_of_reads=$(\
        grep "total reads" "${PROCESSING_DIRECTORY}/model.txt" | \
        grep -Po ": \d+" | \
        grep -Po "\d+" \
    )
}

main() {
    config_file=$1
    validate_config_file "${config_file}"
    source "${config_file}" || exit 1
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
    build_model
    get_fragment_length 
    get_read_length
    get_number_of_reads
    get_coverage_track
    get_bias_track
    if [[ -f "${CONTROL_FILE}" ]]; then
        scale_bias_track
    fi
    get_p_values
    get_cutoff_analysis
    get_best_cutoff
    get_unmerged_peaks
    get_merged_peaks
}

if [[ $#1 -ne 1 ]]; then usage; fi
main "$1"
