#!/bin/bash
#SBATCH --job-name=mh_G1
#SBATCH --partition=general
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=22
#SBATCH --mem=200G
#SBATCH --time=48:00:00
#SBATCH --account=a_zenger
#SBATCH --array=1-20
#SBATCH -o logs/mh_G1_%A_%a.out
#SBATCH -e logs/mh_G1_%A_%a.err

set -euo pipefail

module load anaconda3/2023.09-0
source activate microhaplotype
module load beagle/5.4.22jul22.46e-java-11

cd /scratch/user/aguswibowo/Research/simulation_forward

BASE=/scratch/user/aguswibowo/Research/simulation_forward
OUTDIR=${BASE}/output_G1
NCORES=20
MAP=${BASE}/raw_data/salmon_60k_map_final.txt
REF_VCF=${BASE}/raw_data/salmon_60k_pairs125bp.vcf.gz

ITER=$(printf "%02d" ${SLURM_ARRAY_TASK_ID})
ITER_DIR=${OUTDIR}/iter${ITER}

# ── Runtime + resource CSV (added 2026-05-20) ─────────────────────────────────
# Combined CSV across array tasks via flock-protected append.
# Per step: duration, max RSS, avg RSS, avg CPU%.
# Plus "TOTAL" row per iter at end.
RUNTIME_DIR=${BASE}/runtime
RUNTIME_CSV=${RUNTIME_DIR}/G1_block_def.csv
mkdir -p ${RUNTIME_DIR}

record_step() {
  local step="$1"  start="$2"  end="$3"  status="${4:-0}"
  local max_rss="${5:-0}"  avg_rss="${6:-0}"  avg_cpu="${7:-0}"  n_samples="${8:-0}"
  local duration=$((end - start))
  (
    flock -w 30 200
    if [ ! -s "${RUNTIME_CSV}" ]; then
      echo "iter,step,start_unix,end_unix,duration_sec,exit_status,max_rss_kb,avg_rss_kb,avg_cpu_pct,n_samples" > "${RUNTIME_CSV}"
    fi
    echo "${ITER},${step},${start},${end},${duration},${status},${max_rss},${avg_rss},${avg_cpu},${n_samples}" >> "${RUNTIME_CSV}"
  ) 200>"${RUNTIME_CSV}.lock"
}

# Run a command with /usr/bin/time -v (max RSS + avg CPU%) + background
# RSS poller (avg RSS). Supports pipes via `bash -c "..."` wrapper.
run_timed() {
  local step_name="$1"
  shift
  local time_log=$(mktemp)
  local samples_log=$(mktemp)
  local start=$(date +%s)

  /usr/bin/time -v -o ${time_log} "$@" &
  local cmd_pid=$!

  # Background poller: sample RSS of cmd_pid's child every 10s
  (
    sleep 5
    while kill -0 ${cmd_pid} 2>/dev/null; do
      child=$(pgrep -P ${cmd_pid} 2>/dev/null | head -1)
      target_pid=${child:-${cmd_pid}}
      rss=$(awk '/VmRSS:/ {print $2}' /proc/${target_pid}/status 2>/dev/null)
      [ -n "${rss}" ] && echo "${rss}" >> ${samples_log}
      sleep 10
    done
  ) &
  local poller_pid=$!

  wait ${cmd_pid}
  local exit_code=$?
  wait ${poller_pid} 2>/dev/null

  local end=$(date +%s)

  # Parse /usr/bin/time -v output
  local max_rss=0 avg_cpu=0
  if [ -s "${time_log}" ]; then
    max_rss=$(awk -F': ' '/Maximum resident set size/ {gsub(" ","",$2); print $2}' ${time_log})
    avg_cpu=$(awk -F': ' '/Percent of CPU this job got/ {gsub("[ %]","",$2); print $2}' ${time_log})
  fi

  # Parse poller samples for avg RSS
  local avg_rss=0 n_samples=0
  if [ -s "${samples_log}" ]; then
    read avg_rss n_samples < <(awk '{sum+=$1; n++} END {if (n>0) printf "%d %d\n", sum/n, n; else print "0 0"}' ${samples_log})
  fi

  rm -f ${time_log} ${samples_log}
  record_step "${step_name}" ${start} ${end} ${exit_code} \
              "${max_rss:-0}" "${avg_rss:-0}" "${avg_cpu:-0}" "${n_samples:-0}"
  return ${exit_code}
}

ITER_START=$(date +%s)
echo "Starting MH pipeline iter=${ITER}"

echo "Step 2a: convert-to-vcf iter=${ITER}"
run_timed "convert-to-vcf" ./convert-to-vcf \
    -i ${ITER_DIR}/geno_G1.csv \
    -m ${MAP} \
    -o ${ITER_DIR}/geno_G1.vcf \
    --split-by-chr --parallel --ncores ${NCORES} --merge-after-split

echo "Step 2a.5: Remove duplicate positions iter=${ITER}"
CONTIG_HEADER=${ITER_DIR}/contig_header.txt
bcftools view -h ${REF_VCF} | grep "^##contig" > ${CONTIG_HEADER}

run_timed "dedup" bash -c "bcftools annotate --header-lines ${CONTIG_HEADER} ${ITER_DIR}/geno_G1.vcf | bcftools norm --rm-dup all -O v -o ${ITER_DIR}/geno_G1.vcf.tmp"

mv ${ITER_DIR}/geno_G1.vcf.tmp ${ITER_DIR}/geno_G1.vcf
rm -f ${CONTIG_HEADER}

echo "Step 2b: BEAGLE phasing iter=${ITER}"
run_timed "beagle" java -Xmx50g -jar ${EBROOTBEAGLE}/beagle.jar \
    gt=${ITER_DIR}/geno_G1.vcf \
    nthreads=${NCORES} \
    out=${ITER_DIR}/geno_G1_phased

gunzip -f ${ITER_DIR}/geno_G1_phased.vcf.gz   # -f: overwrite existing .vcf from prior failed attempt

echo "Step 2c: convert-from-vcf iter=${ITER}"
run_timed "convert-from-vcf" ./convert-from-vcf \
    -i ${ITER_DIR}/geno_G1_phased.vcf \
    --map ${MAP} \
    --genofolder ${ITER_DIR}/genoG1 \
    --hapfolder ${ITER_DIR}/haploG1

run_with_retry() {
    local max_attempts=3
    local delay=60
    for attempt in $(seq 1 $max_attempts); do
        "$@" && return 0
        echo "Attempt ${attempt}/${max_attempts} failed. Retrying in ${delay}s..."
        sleep $delay
    done
    echo "All ${max_attempts} attempts failed."
    return 1
}

# Combine retry + monitoring (fix 2026-05-20): /usr/bin/time cannot directly
# execute bash functions, so retry must be the OUTER loop wrapping run_timed
# (which then invokes /usr/bin/time on the actual binary).
run_timed_retry() {
    local step_name="$1"
    shift
    local max_attempts=3 delay=60 attempt
    for attempt in $(seq 1 $max_attempts); do
        run_timed "$step_name" "$@" && return 0
        echo "Attempt ${attempt}/${max_attempts} of '$step_name' failed. Retrying in ${delay}s..."
        sleep $delay
    done
    echo "All ${max_attempts} attempts of '$step_name' failed."
    return 1
}

echo "Step 2d: haplotype-hybrid MH blocks iter=${ITER}"
mkdir -p ${ITER_DIR}/mh_info_ld_haploblock
run_timed_retry "haplotype-hybrid_mh" ./haplotype-hybrid \
    --ncores ${NCORES} \
    -i ${ITER_DIR}/haploG1/chr* \
    -m ${MAP} \
    -o ${ITER_DIR}/mh_info_ld_haploblock \
    --generate-genotypes ${ITER_DIR}/mh_genotypes \
    -v \
    ld-haploblock micro --window-bp 125

echo "Step 2e: haplotype-hybrid ALL blocks iter=${ITER}"
run_timed_retry "haplotype-hybrid_hap" ./haplotype-hybrid \
    --ncores ${NCORES} \
    -i ${ITER_DIR}/haploG1/chr* \
    -m ${MAP} \
    -o ${ITER_DIR}/all_info_haploblock \
    --generate-genotypes ${ITER_DIR}/all_genotypes \
    -v \
    fixed-kb --window-bp 100000

# Per-iter TOTAL summary row (duration only; max/avg RAM aggregated separately
# from per-step rows during post-analysis)
ITER_END=$(date +%s)
record_step "TOTAL" ${ITER_START} ${ITER_END} 0 0 0 0 0

# Validate that G0 QTL block structure is preserved in G1-discovered blocks.
# Reports % QTL blocks intact (all causal SNPs in one G1 block) vs fragmented.
# Output: ${ITER_DIR}/qtl_block_intactness.rds + console summary.
INTACT_CHECK=${BASE}/MicrohapsSel/simulation/post-simulation/check_qtl_block_intactness.R
if [ -f "${INTACT_CHECK}" ]; then
  echo "Step 3: QTL block intactness check iter=${ITER}"
  Rscript "${INTACT_CHECK}" "${SLURM_ARRAY_TASK_ID}" "${OUTDIR}" || \
    echo "  (intactness check failed — pipeline continues, check log)"
fi

echo "Done MH pipeline iter=${ITER}"