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

# Set up files to normalize counts 
CountData <- read.csv("counts_cupcake.mapped_fl_count.txt", row.names = 1)
colnames(CountData) <- c("t0a", "t0c","t10a","t10b","t10c","t2a", "t2b", "t2c","t4a","t4b","t4c")
ColData <- read.table("phenotypes.txt", row.names = 1, header = T)
CountData$total_counts <- rowSums(CountData)

# Bring in sqanti classification file
sqanti <- fread("SQANTI3_results_full_classification.txt")
gencodev38 <- as.data.frame(rtracklayer::import("gencode.v38.annotation.gtf"))

# Keep only genes, remove any other categories (aka exon, transcript etc)
gencodev38 <- filter(gencodev38, type == "gene")
gene_count_lr <- setDT(CountData, keep.rownames = TRUE)[]
gene_count_lr <- inner_join(gene_count_lr, sqanti, by = c("rn" = "isoform"))
gene_count_lr <- left_join(gene_count_lr, gencodev38, by = c("associated_gene" = "gene_id"))

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

#Isoform quantification final criteria (expressed in 25 percent of samples (all replicates of at least on time-point) and has an isoform percentage of 1% or more)
criteria1 <- filter(tpm_normalized_lr, isoin_25 ==1)
criteria1_3 <- filter(criteria1, isoform_percentage >= 0.01) 
criteria1_3 <- criteria1_3[,-c("transcript_name")]

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

# Actual normalization using DESeq2
CountData <- read.csv("counts_cupcake.mapped_fl_count.txt", row.names = 1)
colnames(CountData) <- c("t0a", "t0c","t10a","t10b","t10c","t2a", "t2b", "t2c","t4a","t4b","t4c")
ColData <- read.table("phenotypes.txt", row.names = 1, header = T)
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
