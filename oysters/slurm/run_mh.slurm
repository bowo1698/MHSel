#!/bin/bash
#SBATCH --job-name=masreml_pmx_mh_29m
#SBATCH --qos=normal
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --array=1-7
#SBATCH -o logs/masreml_pmx_mh_29m_%A_%a.out
#SBATCH -e logs/masreml_pmx_mh_29m_%A_%a.err

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export FLEXIBLAS_NUM_THREADS=1
export R_PARALLEL_PORT=random

module load r/4.4.0-gfbf-2023a

echo "========================================"
echo "Job ID: $SLURM_JOB_ID | Array: $SLURM_ARRAY_TASK_ID | Rep: ${REP:-1}"
echo "Node: $(hostname)"
echo "Started: $(date)"
echo "========================================"

BASE="/scratch/user/aguswibowo/Research"
DATA="${BASE}/data/microhaplotypes"
BASE_OUT="${BASE}/prediction_oysters/output4/pmx_mh_29m_iter1"
MODEL="${BASE}/prediction_oysters/model_run"

# REP defaults to 1 when not passed
REP=${REP:-1}
SEED=$((122 + REP))

case $SLURM_ARRAY_TASK_ID in
  1) TRAIT="ShellWeight"; RESPONSE="gaussian"; FIXED="GonCond" ;;
  2) TRAIT="DVH";         RESPONSE="gaussian"; FIXED="GonCond" ;;
  3) TRAIT="APW";         RESPONSE="gaussian"; FIXED="GonCond" ;;
  4) TRAIT="ShellDepth";  RESPONSE="gaussian"; FIXED="GonCond" ;;
  5) TRAIT="colour_Y";    RESPONSE="binary";   FIXED="NULL"    ;;
  6) TRAIT="colour_S";    RESPONSE="binary";   FIXED="NULL"    ;;
  7) TRAIT="colour_Gr";   RESPONSE="binary";   FIXED="NULL"    ;;
esac

if [ "$FIXED" = "NULL" ]; then FIXED_R="NULL"; else FIXED_R="\"${FIXED}\""; fi

if [ "$RESPONSE" = "binary" ]; then
  PI_MH="c(0.995, 0.003, 0.001, 0.001)"
  VC_MH="c(0, 0.0001, 0.001, 0.01)"
else
  PI_MH="c(0.95,0.02,0.02,0.01)"
  VC_MH="c(0, 0.0001, 0.001, 0.01)"
fi

OUT_DIR="${BASE}/prediction_oysters/output4/pmx_mh_29m_iter${REP}/${TRAIT}"

echo "Rep: $REP | Trait: $TRAIT | Seed: $SEED | Out: $OUT_DIR"

TMPCONFIG=$(mktemp /tmp/config_XXXXXX.R)
cat > $TMPCONFIG << EOF
config <- list(
  base_dir          = "${MODEL}",
  data_dir          = "${DATA}",
  output_dir        = "${OUT_DIR}",
  base_output_dir   = "${BASE_OUT}/${TRAIT}",
  block_dir         = "mh_genotypes_filtered_12k",
  phenotype_file    = "phenotype_29m.csv",
  phenotype_trait   = "${TRAIT}",
  fixed_effect      = ${FIXED_R},
  extra_cols        = NULL,
  response_type     = "${RESPONSE}",
  n_chromosomes     = 14,
  min_allele_freq   = 0.01,
  max_allele_freq   = 0.99,
  drop_baseline     = TRUE,
  k_folds           = 5,
  cv_seed           = ${SEED},
  rep_id            = ${REP},
  stratify_by_sire  = FALSE,
  cv_method         = "grm_stratified",
  n_grm_clusters    = 10,
  n_iter            = 50000,
  n_burn            = 25000,
  n_thin            = 10,
  seed              = 123,
  bayesR = list(
    pi        = ${PI_MH},
    var_class = ${VC_MH},
    prior_df  = list(residual = 10, genetic = 10)
  ),
  save_mcmc_samples = FALSE,
  coords_mh_file    = "${DATA}/mh_info_filtered_12k/stats/microhaplotype_coordinates.csv",
  gwas              = list(enabled = TRUE, windsize = 500e3)
)
EOF

echo "--- Trait: $TRAIT | Rep: $REP | Started: $(date) ---"
start_time=$(date +%s)
export GENOMIC_CONFIG=$TMPCONFIG
Rscript ${MODEL}/main_mh.R
echo "Completed in $(($(date +%s) - start_time))s"
rm -f $TMPCONFIG

echo "========================================"
echo "Trait $TRAIT | Rep $REP done: $(date)"
echo "========================================"