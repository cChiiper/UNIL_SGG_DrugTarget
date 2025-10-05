#!/bin/bash -l

# Replace all hard-coded paths (marked with ...)... with actual paths before running.

# List of traits to process – add or change as needed.
declare -a traits=("AD" "AF" "ASTHMA" "BIPOLAR" "BMD" "CAD_consortium" "COPD_FIG_UKBB" "DBP_meta" "ENDOMETRIOSIS_FIG_UKBB" "EPILEPSY" "GLAUCOMA" "IBD_consortium" "IBS" "LDL" "LOAD_consortium" "MDD" "MS" "OA" "PD" "PNEUMONIA_FIG_UKBB" "PSORIASIS_FIG_UKBB" "RA_consortium" "SBP_meta" "SCZ_consortium" "STROKE" "T1D" "T2D" "TC" "VTE_UKBB" "eGFR")

# Input raw GWAS sumstats files are assumed to be here.
SUMSTATS_DIR="/.../"
# Munged files will be stored here (each trait in its own folder).
RESULTS_DIR="/.../"
MERGE_ALLELES="/.../w_hm3.snplist"
LD_CHR_DIR="/.../eur_w_ld_chr/"
FRQ_COL="Freq"
CHUNKSIZE=500000

# Declare an associative array to hold the munge job IDs for each trait.
declare -A munge_job_ids

# Loop over each trait to (re)generate munged sumstats if necessary.
for trait in "${traits[@]}"; do
    input_file="${SUMSTATS_DIR}/${trait}_gwas_summary_uk10kck.ma"
    output_dir="${RESULTS_DIR}/${trait}"
    output_prefix="${output_dir}/${trait}"  # This prefix will be used by munge_sumstats.py.
    munged_file="${output_prefix}.sumstats.gz"  # Adjust extension if needed.
    
    # Create target folder if it does not exist.
    mkdir -p "$output_dir"
    
    if [[ -f "$munged_file" ]]; then
        echo "Munged sumstats for ${trait} already exists. Skipping munge job."
        munge_job_ids[$trait]=""
    else
        echo "Submitting munge job for ${trait}."
        job_id=$(sbatch --parsable <<EOT
#!/bin/bash -l
#SBATCH --job-name=munge_${trait}
#SBATCH --nodes=1
#SBATCH --time=0-00:07:00
#SBATCH --partition=cpu
#SBATCH --mem=16GB
#SBATCH --output=/.../${trait}_%j.out

source /.../bin/activate
conda activate ldsc
/.../ldsc/munge_sumstats.py \
--sumstats ${input_file} \
--out ${output_prefix} \
--merge-alleles ${MERGE_ALLELES} \
--frq ${FRQ_COL} \
--chunksize ${CHUNKSIZE}
EOT
)
        munge_job_ids[$trait]="$job_id"
    fi
done

# Now, loop pairwise over traits to run LDSC genetic correlation.
n=${#traits[@]}
for ((i=0; i<n; i++)); do
    for ((j=i+1; j<n; j++)); do
        trait1=${traits[$i]}
        trait2=${traits[$j]}
        
        file1="${RESULTS_DIR}/${trait1}/${trait1}.sumstats.gz"
        file2="${RESULTS_DIR}/${trait2}/${trait2}.sumstats.gz"
        output_prefix_pair="${RESULTS_DIR}/${trait1}_${trait2}"
        
        # Build dependency option if any munge jobs were submitted.
        dependency=""
        dep_ids=()
        if [[ -n "${munge_job_ids[$trait1]}" ]]; then
            dep_ids+=("${munge_job_ids[$trait1]}")
        fi
        if [[ -n "${munge_job_ids[$trait2]}" ]]; then
            dep_ids+=("${munge_job_ids[$trait2]}")
        fi
        if [ ${#dep_ids[@]} -gt 0 ]; then
            # Join job IDs with ':' for SLURM dependency.
            dependency="--dependency=afterok:$(IFS=:; echo "${dep_ids[*]}")"
        fi
        
        sbatch $dependency <<EOT
#!/bin/bash -l
#SBATCH --job-name=LDSC_${trait1}_${trait2}
#SBATCH --nodes=1
#SBATCH --time=0-00:05:00
#SBATCH --partition=cpu
#SBATCH --mem=4GB
#SBATCH --output=/.../LDSC_${trait1}_${trait2}_%j.out

source /.../bin/activate
conda activate ldsc
cd /.../

/.../ldsc/ldsc.py \
--rg ${file1},${file2} \
--ref-ld-chr ${LD_CHR_DIR}/ \
--w-ld-chr ${LD_CHR_DIR}/ \
--out ${output_prefix_pair}
EOT

    done
done
