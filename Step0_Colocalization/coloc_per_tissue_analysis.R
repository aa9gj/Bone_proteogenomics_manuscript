# Capture Command Line Arguments
cmdArgs <- commandArgs(trailingOnly = TRUE)
morris <- cmdArgs[2]
gtex <- cmdArgs[3]
snps <- cmdArgs[4]
morris_id <- cmdArgs[5]
gtex_id <- cmdArgs[6]
coloc_res <- cmdArgs[7]

# Using sQTL colocalization with coloc
library(data.table)
library(dplyr)
library(tidyr)
library(tidyverse)
library(coloc)

morris_id <- fread(file = morris, header = FALSE)
morris_id <- morris_id[!duplicated(morris_id$V7),]

gtex_id <- fread(file = gtex, header = FALSE)

setwd("/scratch/aa9gj/wrk/sQTL_GTEx_Data/GTEx_sqtl_data/analysis_output_by_tissue/snps/Cells_Cultured_fibroblasts_snps")
my_files <- fread(snps,header=FALSE)
my_data <- list()
for (i in seq_along(my_files$V1)) {
	my_data[[i]] <- fread(file = my_files$V1[i], header =FALSE)
	my_data[[i]] <- subset(my_data[[i]], my_data[[i]]$V3 < 200000 & my_data[[i]]$V3 > -200000)
	my_data[[i]] <- inner_join(my_data[[i]], gtex_id, by = c("V2"="V1"))
	my_data[[i]] <- my_data[[i]][!duplicated(my_data[[i]]$V7.y),]
	my_data[[i]] <- subset(my_data[[i]], my_data[[i]]$V6.x > 0)
}

my_coloc <- list()
for (i in seq_along(my_files$V1)) {
 my_coloc[[i]] <- coloc.abf(dataset1=list(snp=my_data[[i]]$V7.y, type ="quant", MAF=my_data[[i]]$V6.x, N=500, pvalues = my_data[[i]]$V7.x), dataset2=list(snp= morris_id$V7, type ="quant", MAF= morris_id$V18, N=426824, pvalues = morris_id$V16))
}

# Extracting necessary coloc results
summary_list <- list()
for (i in seq_along(my_files$V1)) {
 summary_list[[i]] <- my_coloc[[i]]$summary
}

df <- as.data.frame(do.call("rbind", summary_list))
df$event <- my_files$V1
setwd("/scratch/aa9gj/wrk/sQTL_GTEx_Data/GTEx_sqtl_data/analysis_output_by_tissue/coloc_analysis/Cells_Cultured_fibroblasts_coloc")
write.table(df, file = coloc_res, row.names = FALSE, col.names = FALSE, quote=FALSE, sep = "\t")
