#!/bin/bash

#SBATCH -J tissue
#SBATCH -n 1
#SBATCH -t 72:00:00
#SBATCH -o coloc.out
#SBATCH -e coloc.err
#SBATCH -p standard
#SBATCH -A cphg-farber
#SBATCH --mem=100000
#SBATCH --array=1-22

i=$SLURM_ARRAY_TASK_ID
line=$(sed -n "${i}{p;}" < morris_snp_paths) #Where eBMD SNPs files are by chromosome 
mine=$(sed -n "${i}{p;}" < gtex_files_path) #Where GTEx files are located
fine=$(sed -n "${i}{p;}" < snps/Cells_Cultured_fibroblasts_snps/snp_files) #list of snp file names
wine=$(sed -n "${i}{p;}" < chr_ids) #list of chromosme ids (chr1 etc)
sine=$(sed -n "${i}{p;}" < current_gtex_ids) #list of GTEx chr ids (chr1 etc)
kine=$(sed -n "${i}{p;}" < coloc_list) #list coloc chr ids (chr1 etc)

module purge
module load goolf/7.1.0_3.1.4 R/4.0.3
Rscript coloc_per_tissue_analysis.R ${SLURM_ARRAY_TASK_ID} ${line} ${mine} ${fine} ${wine} ${sine} ${kine}
