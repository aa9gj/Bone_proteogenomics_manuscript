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
line=$(sed -n "${i}{p;}" < morris_snp_paths)
mine=$(sed -n "${i}{p;}" < gtex_files_path)
fine=$(sed -n "${i}{p;}" < /scratch/aa9gj/wrk/sQTL_GTEx_Data/GTEx_sqtl_data/analysis_output_by_tissue/snps/Cells_Cultured_fibroblasts_snps/snp_files)
wine=$(sed -n "${i}{p;}" < chr_ids)
sine=$(sed -n "${i}{p;}" < current_gtex_ids)
kine=$(sed -n "${i}{p;}" < coloc_list)

module purge
module load goolf/7.1.0_3.1.4 R/4.0.3
Rscript coloc_per_tissue_analysis.R ${SLURM_ARRAY_TASK_ID} ${line} ${mine} ${fine} ${wine} ${sine} ${kine}
