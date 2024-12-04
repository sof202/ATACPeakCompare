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
    get_fragment_length 
    get_read_length
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
