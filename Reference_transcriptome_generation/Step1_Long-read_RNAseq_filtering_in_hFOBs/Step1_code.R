###########################################################################
###########################################################################
###                                                                     ###
###        Long-read RNAseq filtering in hFOBs                          ###     
### 								                ###
###########################################################################
###########################################################################

# Where do I get these files to run this code? 
# Download the following files from the Zenodo repo (counts_cupcake.mapped_fl_count.txt; phenotypes.txt; SQANTI3_results_full_classification.txt)
# Download GENCODE v38 gtf file from GENCODE

# Set working directory to script location
script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(script_dir)

library(dplyr)
library(tibble)
library(data.table)
library(stringr)

# Set up files to normalize counts 
CountData <- read.csv("counts_cupcake.mapped_fl_count.txt", row.names = 1)
colnames(CountData) <- c("t0a", "t0c","t10a","t10b","t10c","t2a", "t2b", "t2c","t4a","t4b","t4c")
ColData <- read.csv("phenotype.csv", row.names = 1, header = T)
CountData$total_counts <- rowSums(CountData)

# Bring in sqanti classification file
sqanti <- fread("SQANTI3_results_full_classification.txt")
gencodev38 <- as.data.frame(rtracklayer::import("gencode.v38.annotation.gtf"))

# Keep only genes, remove any other categories (aka exon, transcript etc)
gencodev38 <- filter(gencodev38, type == "gene")
gene_count_lr <- setDT(CountData, keep.rownames = TRUE)[]
gene_count_lr <- inner_join(gene_count_lr, sqanti, by = c("rn" = "isoform"))
gene_count_lr <- left_join(gene_count_lr, gencodev38, by = c("associated_gene" = "gene_id"))


# GS TRBL
# countdata_only_pbs <- setdiff(CountData$rn, sqanti$isoform)
# CountData %>% filter(rn %in% countdata_only_pbs)
# setdiff(sqanti$isoform, CountData$rn)

# Remove novel genes which are most likely artifacts
gene_count_lr <- gene_count_lr[!is.na(gene_count_lr$gene_name),]
gene_count_lr <- filter(gene_count_lr, gene_type == "protein_coding")

# Doing TPM manually for total gene counts (for filteration purposes)
gene_count_lr$step1 <- gene_count_lr$total_counts / gene_count_lr$length
gene_count_lr_matrix <- gene_count_lr[,c(1,85)]
gene_count_lr_matrix <- gene_count_lr_matrix %>% remove_rownames %>% column_to_rownames(var="rn")
gene_count_lr_matrix <- as.data.frame(t( t(gene_count_lr_matrix) * 1e6 / colSums(gene_count_lr_matrix) ))
tpm_normalized_lr <- setDT(gene_count_lr_matrix, keep.rownames = TRUE)[]
colnames(tpm_normalized_lr) <- c("rn", "total_tpm")
tpm_normalized_lr <- inner_join(gene_count_lr, tpm_normalized_lr, by = "rn")

# Create isoforms percentages to use as a threshold with TPM 
isoform_percent <- aggregate(tpm_normalized_lr$total_tpm, by=list(Category=tpm_normalized_lr$gene_name), FUN=sum)
tpm_normalized_lr <- left_join(tpm_normalized_lr, isoform_percent, by = c("gene_name" = "Category"))
tpm_normalized_lr$isoform_percentage <- tpm_normalized_lr$total_tpm/tpm_normalized_lr$x

# Add a column of to ensure that the isoform is expressed in all replicates of a time point, >0 meaning true, Inf meaning false
tpm_normalized_lr$t0a_exp <- 1/tpm_normalized_lr$t0a
tpm_normalized_lr$t0c_exp <- 1/tpm_normalized_lr$t0c
tpm_normalized_lr$t2a_exp <- 1/tpm_normalized_lr$t2a
tpm_normalized_lr$t2b_exp <- 1/tpm_normalized_lr$t2b
tpm_normalized_lr$t2c_exp <- 1/tpm_normalized_lr$t2c
tpm_normalized_lr$t4a_exp <- 1/tpm_normalized_lr$t4a
tpm_normalized_lr$t4b_exp <- 1/tpm_normalized_lr$t4b
tpm_normalized_lr$t4c_exp <- 1/tpm_normalized_lr$t4c
tpm_normalized_lr$t10a_exp <- 1/tpm_normalized_lr$t10a
tpm_normalized_lr$t10b_exp <- 1/tpm_normalized_lr$t10b
tpm_normalized_lr$t10c_exp <- 1/tpm_normalized_lr$t10c

# Create a column to determine if the isoform is expressed in 25% of samples (all replicates of at least on time-point)
# 2023-12-23 GS added this as guess because code missing from AA original.
tpm_normalized_lr <- tpm_normalized_lr %>% mutate(isoin_25 = as.integer((t0a > 0 & t0c > 0) | (t2a > 0 & t2b > 0 & t2c > 0) | (t4a > 0 & t4b > 0 & t4c > 0) | (t10a > 0 & t10b > 0 & t10c > 0)))

#Isoform quantification final criteria (expressed in 25 percent of samples (all replicates of at least on time-point) and has an isoform percentage of 1% or more)
criteria1 <- filter(tpm_normalized_lr, isoin_25 ==1)
criteria1_3 <- filter(criteria1, isoform_percentage >= 0.01) 
criteria1_3 <- criteria1_3[,-c("transcript_name")]

# GS TRBL - discrepancy in final table
# # Merge the counts
# merged_counts <- full_join(aa_counts, criteria_counts, by = "gene_name")
# 
# # Identify gene_names where counts differ
# genes_with_diff_counts <- merged_counts %>%
#   filter(is.na(aa_count) | is.na(criteria_count) | aa_count != criteria_count) %>%
#   pull(gene_name)
# 
# # Extract rows from both data frames for these genes
# aa_rows_diff <- aa_tpm_norm_table %>% filter(gene_name %in% genes_with_diff_counts) %>% select(rn, gene_name) %>% mutate(source="aa") 
# criteria_rows_diff <- criteria1_3 %>% filter(gene_name %in% genes_with_diff_counts) %>% select(rn, gene_name) %>% mutate(source="criteria")
# 
# # Combine the results (optional)
# combined_diff <- bind_rows(aa_rows_diff, criteria_rows_diff)
# 
# 
# 
# col_names <- grep("^t\\d+[a-c]$", colnames(aa_tpm_norm_table), value=TRUE)
# mismatched <- anti_join(aa_tpm_norm_table, criteria1_3, by = col_names)
# 
# unique_to_aa <- anti_join(aa_tpm_norm_table, criteria1_3, by = "rn")
# unique_to_criteria1_3 <- anti_join(criteria1_3, aa_tpm_norm_table, by = "rn")
# unique_to_aa <-  unique_to_aa %>% mutate(source="aa_table") %>% select(source, everything())
# unique_to_criteria1_3 <- unique_to_criteria1_3 %>% mutate(source="criteria_table") %>% select(source, everything())
# mismatched_rn <- bind_rows(unique_to_aa, unique_to_criteria1_3)
# 
# unique_genes_in_mismatch <- unique(c(unique_to_aa$gene_name, unique_to_criteria1_3$gene_name))
# 
# output_file_path <- "output.txt"
# for (gene in unique_genes_in_mismatch) {
#   cat("*****\n", file = output_file_path, append = TRUE)
#   write.table(aa_tpm_norm_table %>% filter(gene_name == gene) %>% select(gene_name, 1:5, "length"), file = output_file_path, append = TRUE, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
#   cat("*\n", file = output_file_path, append = TRUE)
#   write.table(criteria1_3 %>% filter(gene_name == gene) %>% select(gene_name, 1:5, "length"), file = output_file_path, append = TRUE, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
#   cat("*****\n", file = output_file_path, append = TRUE)
# }
# 
# aa_tpm_norm_table %>% filter(gene_name==gene) %>% select("gene_name", 1:5, "length")
# criteria1_3 %>% filter(gene_name==gene) %>% select("gene_name", 1:5, "length")
# 
# table(criteria1_3$structural_category)

# Write out table of Criteria 1/3 filtered hfob isoforms (reported in main manuscript)


# Summary data for long read isoforms (if interested)
median(criteria1_3$length)
min(criteria1_3$length)
max(criteria1_3$length)
range(criteria1_3$length)
novelty <- filter(criteria1_3, associated_transcript == "novel")
mean(novelty$length)
median(novelty$length)
not_novel <- filter(criteria1_3, associated_transcript != "novel")
mean(not_novel$length)
median(not_novel$length)

# Counts for structural categories (These structural categories are described fully in the SQANTI paper)
FSM <- filter(criteria1_3, structural_category == "full-splice_match")
ISM <- filter(criteria1_3, structural_category == "incomplete-splice_match")
genic <- filter(criteria1_3, structural_category == "genic")
NIC <- filter(criteria1_3, structural_category == "novel_in_catalog")
NNC <- filter(criteria1_3, structural_category == "novel_not_in_catalog")

# Number known isoforms
num_known <- nrow(FSM) + nrow(ISM)

# Number of novel isoforms
num_novel <- nrow(NIC) + nrow(NNC)

# Fraction novel
num_novel/(num_known + num_novel)


# Actual normalization using DESeq2
CountData <- read.csv("counts_cupcake.mapped_fl_count.txt", row.names = 1)
colnames(CountData) <- c("t0a", "t0c","t10a","t10b","t10c","t2a", "t2b", "t2c","t4a","t4b","t4c")
ColData <- read.csv("phenotype.csv", row.names = 1, header = T)
reorder_idx <- match(rownames(ColData),colnames(CountData))
# Reorder the columns of the count data
CountData <- CountData[ , reorder_idx]
all(rownames(ColData) %in% colnames(CountData))
cpm <- DESeqDataSetFromMatrix(countData =  CountData,
                              colData =  ColData,
                              design = ~ Condition)
cpm <- estimateSizeFactors(cpm)

# Extract the normalized counts
cpm_normalized <- as.data.frame(counts(cpm, normalize = TRUE))
cpm_normalized <- tibble::rownames_to_column(cpm_normalized, "rn")
colnames(cpm_normalized) <- c("rn","t0a_norm", "t0c_norm", "t2a_norm", "t2b_norm", "t2c_norm", "t4a_norm", "t4b_norm", "t4c_norm", "t10a_norm", "t10b_norm", "t10c_norm")

# Add normalized counts information to the thresholded file
criteria1_3 <- inner_join(cpm_normalized, criteria1_3, by = "rn")

# We end up with this number of genes
length(unique(criteria1_3$gene_name))#12068 protein-coding genes 

