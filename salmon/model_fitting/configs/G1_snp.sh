#!/bin/bash
#SBATCH --job-name=fit_cv_G1_snp
#SBATCH --array=1-20 #1-20
#SBATCH --partition=general
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=200G
#SBATCH --time=120:00:00
#SBATCH --account=a_zenger
#SBATCH -o logs/fit_cv_G1_snp_%A_%a.out
#SBATCH -e logs/fit_cv_G1_snp_%A_%a.err

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export FLEXIBLAS_NUM_THREADS=1

module load r/4.4.0-gfbf-2023a

ITER=$SLURM_ARRAY_TASK_ID
ITERPAD=$(printf '%02d' $ITER)
SEED=$((ITER * 42))
MODEL="/scratch/user/aguswibowo/Research/predict_salmon_forward/model_run/test7"
DATA_DIR_G1="/scratch/user/aguswibowo/Research/simulation_forward2/output_G1/iter${ITERPAD}"
OUT_BASE="/scratch/user/aguswibowo/Research/predict_salmon_forward/output2/predict_within_gen_G1_full_test7/iter${ITERPAD}"

SCENARIOS_G1='list(
    S1_equal   = "pheno_S1_G1.csv",
    S2_unequal = "pheno_S2_G1.csv",
    S3_equal   = "pheno_S3_G1.csv",
    S4_unequal = "pheno_S4_G1.csv",
    S5_equal   = "pheno_S5_G1.csv",
    S6_unequal = "pheno_S6_G1.csv"
  )'

echo "=== ITER=${ITER} SNP Started: $(date) ==="

TMPCONFIG=$(mktemp /tmp/config_snp_cv_XXXXXX.R)
cat > $TMPCONFIG << EOF
config <- list(
  base_dir            = "${MODEL}",
  data_dir            = "${DATA_DIR_G1}",
  output_dir          = "${OUT_BASE}/snp",
  genotype            = "geno_G1.csv",
  pedigree            = "pedigree_G1.csv",
  population_type     = "reference",
  min_allele_freq     = 0.01,
  callrate_threshold  = 0.8,
  seed                = 123,
  marker_type         = "snp",
  snp_encoding        = "zscore",
  n_iter              = 50000,
  n_burn              = 25000,
  n_thin              = 10,
  bayesR = list(
    pi        = c(0.90, 0.05, 0.03, 0.02),
    var_class = c(0, 0.0001, 0.001, 0.01),
    prior_df  = list(residual = -1, genetic = 2)
  ),
  save_mcmc_samples   = FALSE,
  mcmc_diag           = list(enabled = TRUE, k_window = 20L),
  scenarios           = ${SCENARIOS_G1},
  qtl_dir             = "/scratch/user/aguswibowo/Research/simulation_forward2/output_G1/rds",
  snp_map_file        = "/scratch/user/aguswibowo/Research/simulation_forward2/raw_data/salmon_60k_map_final_v4.txt",
  mh_coords_file      = "/scratch/user/aguswibowo/Research/simulation_forward2/raw_data/mh_info_ld_haploblock_G0/stats/microhaplotype_coordinates.csv",
  hap_coords_file     = "/scratch/user/aguswibowo/Research/simulation_forward2/raw_data/all_info_haploblock_G0/stats/microhaplotype_coordinates.csv",
  gwas                = list(enabled = TRUE, windsize = 500e3)
)
EOF
export GENOMIC_CONFIG=$TMPCONFIG
Rscript ${MODEL}/main-snp.R
rm -f $TMPCONFIG

echo "=== ITER=${ITER} SNP Completed: $(date) ==="