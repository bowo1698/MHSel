#!/bin/bash
#SBATCH --job-name=G0_create_map
#SBATCH --partition=general
#SBATCH --qos=short
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --account=a_zenger
#SBATCH -o logs/G0_create_map_%j.out
#SBATCH -e logs/G0_create_map_%j.err

module load anaconda3/2023.09-0
source activate microhaplotype

DIR=/scratch/user/aguswibowo/Research/simulation_forward/raw_data
VCF=$DIR/salmon_60k_pairs125bp.vcf.gz

bcftools query -f '%ID\t%CHROM\t%POS\n' $VCF > $DIR/salmon_60k_map.txt
echo -e "SNPID\tChr\tPosition" | cat - $DIR/salmon_60k_map.txt > $DIR/salmon_60k_map_final.txt

bgzip -d -k $VCF