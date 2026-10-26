#!/bin/bash
#SBATCH --job-name=masreml_pmx_snp_29m
#SBATCH --qos=normal
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --account=a_jcu_hpc
#SBATCH --array=1-28
#SBATCH -o logs/masreml_pmx_snp_29m_%A_%a.out
#SBATCH -e logs/masreml_pmx_snp_29m_%A_%a.err

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
BASE_OUT="${BASE}/prediction_oysters/output4/pmx_snp_29m_iter1"
MODEL="${BASE}/prediction_oysters/model_run"

REP=${REP:-1}
SEED=$((122 + REP))

PANELS=("snp_target" "snp_full" "snp_rand5k" "snp_from_mh")
GENO_FILES=(
  "snp_target_genotypes.csv"
  "snp_all_filtered_genotypes.csv"
  "snp_all_filtered_rand_5k_genotypes.csv"
  "genotypes_microhap_12k.csv"
)
# snp_target has its own ID namespace; the other three panels are subsets
# of the 12K MH-derived SNP set.
SNP_MAPS=(
  "target_snps_map.txt"
  "map_final_12K.txt"
  "map_final_12K.txt"
  "map_final_12K.txt"
)
TRAITS=("ShellWeight" "DVH" "APW" "ShellDepth" "colour_Y" "colour_S" "colour_Gr")
RESPONSES=("gaussian" "gaussian" "gaussian" "gaussian" "binary" "binary" "binary")
FIXED_EFFECTS=("GonCond" "GonCond" "GonCond" "GonCond" "NULL" "NULL" "NULL")

PANEL_IDX=$(( (SLURM_ARRAY_TASK_ID - 1) / 7 ))
TRAIT_IDX=$(( (SLURM_ARRAY_TASK_ID - 1) % 7 ))

PANEL=${PANELS[$PANEL_IDX]}
GENO_FILE=${GENO_FILES[$PANEL_IDX]}
SNP_MAP_FILE=${SNP_MAPS[$PANEL_IDX]}
TRAIT=${TRAITS[$TRAIT_IDX]}
RESPONSE=${RESPONSES[$TRAIT_IDX]}
FIXED=${FIXED_EFFECTS[$TRAIT_IDX]}

if [ "$FIXED" = "NULL" ]; then FIXED_R="NULL"; else FIXED_R="\"${FIXED}\""; fi

if [ "$RESPONSE" = "binary" ]; then
  PI="c(0.995,0.003,0.001,0.001)"
  VC="c(0, 0.0001, 0.001, 0.01)"
else
  PI="c(0.95,0.02,0.02,0.01)"
  VC="c(0, 0.001, 0.01, 0.1)"
fi

OUT_DIR="${BASE}/prediction_oysters/output4/pmx_snp_29m_iter${REP}/${PANEL}/${TRAIT}"

echo "Rep: $REP | Panel: $PANEL | Trait: $TRAIT | Seed: $SEED"

TMPCONFIG=$(mktemp /tmp/config_XXXXXX.R)
cat > $TMPCONFIG << EOF
config <- list(
  base_dir          = "${MODEL}",
  data_dir          = "${DATA}",
  output_dir        = "${OUT_DIR}",
  base_output_dir   = "${BASE_OUT}/${PANEL}/${TRAIT}",
  genotype_file     = "${GENO_FILE}",
  phenotype_file    = "phenotype_29m.csv",
  phenotype_trait   = "${TRAIT}",
  fixed_effect      = ${FIXED_R},
  extra_cols        = NULL,
  response_type     = "${RESPONSE}",
  min_allele_freq   = 0.01,
  max_allele_freq   = 0.99,
  callrate_threshold = 0.95,
  k_folds           = 5,
  cv_seed           = ${SEED},
  rep_id            = ${REP},
  cv_method         = "grm_stratified",
  n_grm_clusters    = 10,
  marker_type       = "snp",
  snp_encoding      = "zscore",
  n_iter            = 50000,
  n_burn            = 25000,
  n_thin            = 10,
  seed              = 123,
  bayesR = list(
    pi        = ${PI},
    var_class = ${VC},
    prior_df  = list(residual = -1, genetic = 2)
  ),
  save_mcmc_samples = FALSE,
  snp_map_file      = "${DATA}/${SNP_MAP_FILE}",
  gwas              = list(enabled = TRUE, windsize = 500e3)
)
EOF

echo "--- Panel: $PANEL | Trait: $TRAIT | Rep: $REP | Started: $(date) ---"
start_time=$(date +%s)
export GENOMIC_CONFIG=$TMPCONFIG
Rscript ${MODEL}/main_snp.R
echo "Completed in $(($(date +%s) - start_time))s"
rm -f $TMPCONFIG

echo "========================================"
echo "Panel $PANEL | Trait $TRAIT | Rep $REP done: $(date)"
echo "========================================"