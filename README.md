# MHSel

Genomic prediction using **microhaplotypes (MH)** and **SNPs**, run on two datasets:

- `salmon/` — simulated data (AlphaSimR). Forward simulation of G1 from a G0 SNP panel, then prediction.
- `oysters/` — real data. Cross-validated prediction of shell and colour traits.

Each dataset fits the same two models — GBLUP and BayesR — to both marker types, so MH and SNP accuracies are directly comparable.

## Layout

```
salmon/
  pipeline/                 G0/G1 haploblock definition, map creation (bash)
  simulation/simulationG1.R AlphaSimR forward simulation
  model_fitting/            entry points, R/, configs/ (SLURM + config)
oysters/
  main_mh.R, main_snp.R     entry points
  R/                        data prep, W matrix, CV setup, models/
  slurm/                    SLURM jobs (one array task per trait)
```

## Usage

Every entry script reads its config from the `GENOMIC_CONFIG` environment variable
(defaults to `config.R`). The SLURM scripts generate that config into a temp file
and export it, so in practice you submit a job rather than call R directly:

```bash
sbatch salmon/model_fitting/configs/G1_mh.sh         # salmon G1, MH
sbatch oysters/slurm/run_mh.slurm                    # oysters, MH
```

To run outside SLURM, write a `config` list matching the one in the `.slurm` /
`.sh` file, then:

```bash
GENOMIC_CONFIG=/path/to/config.R Rscript oysters/main_mh.R
```

Paths in the configs point at `/scratch/...` on the HPC — edit them for a
different machine.

## Dependencies

R ≥ 4.4. Core:

- [`maspipeline`](https://github.com/bowo1698/maspipeline) — Rust CLI for phasing
  and microhaplotype discovery. `salmon/pipeline/` calls its three binaries
  (`convert-to-vcf`, `convert-from-vcf`, `haplotype-hybrid`) from the working
  directory, so build it and put them on `PATH` or symlink them there.
- [`masbayes`](https://github.com/bowo1698/masbayes) — BayesR / Bayesian marker models
- [`masreml`](https://github.com/bowo1698/masreml) — REML variance components, GBLUP

Also: `tidyverse`, `coda`, `sommer`, `pROC`, `foreach`, `doParallel`,
`RhpcBLASctl`, `data.table`; plus `AlphaSimR` and `vcfR` for the salmon
simulation.

External tools for `salmon/pipeline/`: `bcftools`, `bgzip`, and BEAGLE 5.4
(for G1 phasing).
