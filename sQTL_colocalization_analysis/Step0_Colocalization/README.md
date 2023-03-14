### Step 0: Generation of disease relevant sQTL associations

1. Download sQTL V8 data for all associations from GTEx using https://console.cloud.google.com/storage/browser/gtex-resources

2. Download latest eBMD GWAS summary stats from GEFOS using http://www.gefos.org/?q=content/data-release-2018

3. The eBMD snp coordinates are in hg19, therefore, LiftOver to hg38 using LiftOver in R. 

4. For each lead GWAS association, we identified that the genes (and their sQTL junctions) that are within 400Kb window (200Kb upstream and 200Kb downstream)

5. Then for each sQTL junction, we take all the snps that are within 400Kb window

6. Perform colocalization analysis

### To make things run smoothly we suggest the following: 

1. Slice eBMD GWAS data by chromosome, GTEx look up table by chromosome, and GWAS lead snps by chromosome for computation efficiency. This data can be found already processed on zenodo

2. Create separate directories for each tissue, and make a list of those names. Also, make a list of GTEx file names, and run automate_file_prod_chr[1-22].sh 

3. To filter only genes (and their junctions) that are within 400Kb of GWAS loci using test_auto_gene_prod.R. HOWEVER we need to run make_gene_file_list_test.slurm to produce a prerequesite file.

4. To get the snps that are within a 400kb window of each gene efficiently, we created this script for each tissue SNP_retrieval_v2.py. You can automate it as a slurm to run it for all tissues at the same time. 

5. We recommend separating snps files based on chromosome to be able to use coloc to perform the colocalization analysis using coloc_per_tissue_analysis.R wrapped in coloc_per_tissue_analysis.slurm

