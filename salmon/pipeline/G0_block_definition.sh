#!/bin/bash
#SBATCH --job-name=G0_block_def
#SBATCH --partition=general
#SBATCH --qos=short
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --account=a_zenger
#SBATCH -o logs/G0_block_def_%j.out
#SBATCH -e logs/G0_block_def_%j.err

BASE=/scratch/user/aguswibowo/Research/simulation_forward
RAW=$BASE/raw_data
MAP=$RAW/salmon_60k_map_final.txt
NCORES=15

echo "convert-from-vcf G0"
./convert-from-vcf \
    -i $RAW/salmon_60k_pairs125bp.vcf \
    --map $MAP \
    --genofolder $RAW/genoG0 \
    --hapfolder $RAW/haploG0

echo "haplotype-hybrid MH blocks"
./haplotype-hybrid \
    --ncores $NCORES \
    -i $RAW/haploG0/chr* \
    -m $MAP \
    -o $RAW/mh_info_ld_haploblock_G0 \
    --generate-genotypes $RAW/mh_genotypes_G0 \
    -v \
    ld-haploblock micro --window-bp 125

echo "haplotype-hybrid haplotype blocks"
./haplotype-hybrid \
    --ncores $NCORES \
    -i $RAW/haploG0/chr* \
    -m $MAP \
    -o $RAW/all_info_haploblock_G0 \
    --generate-genotypes $RAW/all_genotypes_G0 \
    -v \
    fixed-kb --window-bp 100000

echo "Done MH pipeline G0"