#!/bin/bash
#SBATCH --job-name=mh_G1
#SBATCH --partition=general
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=22
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --account=a_jcu_hpc
#SBATCH --array=1-20
#SBATCH -o logs/mh_G1_%A_%a.out
#SBATCH -e logs/mh_G1_%A_%a.err

set -euo pipefail

# Environment (modules, RUN_ROOT, RAW_DIR, PATH with maspipeline bins) comes
# from hpc/env.sh, sourced by hpc/salmon_02_blocks.slurm.
: "${RUN_ROOT:?source hpc/env.sh first}" "${RAW_DIR:?source hpc/env.sh first}"

BASE=${RUN_ROOT}/salmon
OUTDIR=${BASE}/output_G1
NCORES=${NCORES:-20}
MAP=${RAW_DIR}/salmon_60k_map_final_v4.txt
REF_VCF=${RAW_DIR}/salmon_60k_pairs125bp.vcf.gz

ITER=$(printf "%02d" "${1:-${SLURM_ARRAY_TASK_ID}}")
ITER_DIR=${OUTDIR}/iter${ITER}

# ── Runtime + resource CSV ────────────────────────────────────────────────────
# One CSV shared across array tasks, appended under flock.
# Per step: duration, max RSS, avg RSS, avg CPU%; plus a TOTAL row per iter.
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

# Run a command under /usr/bin/time -v (max RSS, avg CPU%) with a background
# RSS poller for avg RSS. Pipes need a `bash -c "..."` wrapper.
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
run_timed "convert-to-vcf" convert-to-vcf \
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
run_timed "beagle" java -Xmx50g -jar ${BEAGLE_JAR:-${EBROOTBEAGLE}/beagle.jar} \
    gt=${ITER_DIR}/geno_G1.vcf \
    nthreads=${NCORES} \
    out=${ITER_DIR}/geno_G1_phased

gunzip -f ${ITER_DIR}/geno_G1_phased.vcf.gz   # -f: overwrite .vcf from a failed run

echo "Step 2c: convert-from-vcf iter=${ITER}"
run_timed "convert-from-vcf" convert-from-vcf \
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

# /usr/bin/time cannot execute bash functions, so retry has to be the OUTER
# loop wrapping run_timed, which then times the actual binary.
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
run_timed_retry "haplotype-hybrid_mh" haplotype-hybrid \
    --ncores ${NCORES} \
    -i ${ITER_DIR}/haploG1/chr* \
    -m ${MAP} \
    -o ${ITER_DIR}/mh_info_ld_haploblock \
    --generate-genotypes ${ITER_DIR}/mh_genotypes \
    -v \
    ld-haploblock micro --window-bp 125

echo "Step 2e: haplotype-hybrid ALL blocks iter=${ITER}"
run_timed_retry "haplotype-hybrid_hap" haplotype-hybrid \
    --ncores ${NCORES} \
    -i ${ITER_DIR}/haploG1/chr* \
    -m ${MAP} \
    -o ${ITER_DIR}/all_info_haploblock \
    --generate-genotypes ${ITER_DIR}/all_genotypes \
    -v \
    fixed-kb --window-bp 100000

# Per-iter TOTAL row: duration only; RAM is aggregated from per-step rows later
ITER_END=$(date +%s)
record_step "TOTAL" ${ITER_START} ${ITER_END} 0 0 0 0 0

# Check G0 QTL block structure survives in G1-discovered blocks: % intact
# (all causal SNPs in one G1 block) vs fragmented.
# Writes ${ITER_DIR}/qtl_block_intactness.rds.
INTACT_CHECK=${CODE}/salmon/tools/check_qtl_block_intactness.R   # optional; skipped if absent
if [ -f "${INTACT_CHECK}" ]; then
  echo "Step 3: QTL block intactness check iter=${ITER}"
  Rscript "${INTACT_CHECK}" "${SLURM_ARRAY_TASK_ID}" "${OUTDIR}" || \
    echo "  (intactness check failed — pipeline continues, check log)"
fi

echo "Done MH pipeline iter=${ITER}"