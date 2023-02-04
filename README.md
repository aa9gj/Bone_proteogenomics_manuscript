# Bone_proteogenomics_manuscript
Code accompanying the manuscript

Step0_notes
1) Download sQTL V8 data for all associations from GTEx using https://console.cloud.google.com/storage/browser/gtex-resources 

2) Download latest eBMD GWAS summary stats from GEFOS using http://www.gefos.org/?q=content/data-release-2018

3) Slice eBMD GWAS data by chromosome, GTEx look up table by chromosome, and GWAS lead snps by chromosome for computation efficiency 

4) Create separate directories for each tissue, and make a list of those names. Also, make a list of GTEx file names, and run automate_file_prod_chr[1-22].sh

5) This should produce all the genes that are within eBMD GWAS loci, next we need to filter only genes that are within 200Kb of GWAS loci using test_auto_gene_prod.R wrapped in filter_genes_generated_automated.slurm. HOWEVER we need to run make_gene_file_list_test.slurm to produce a prerequesite file.

6) Now you have the genes that are within a 400kb window of each eBMD GWAS locus. We have to get the snps that are within a 400kb window of each gene using this script for each tissue test_gene_snps.sh 

7) Once that's done, separate those snps files based on chromosome to be able to use coloc to perform the colocalization analysis using test_auto_coloc.R wrapped in run_auto_coloc.slurm

8) Optional step: make sense of the events by visualizaing via a mirrorplot or relating them to gene names using mirrorplot_annotate_events.R
