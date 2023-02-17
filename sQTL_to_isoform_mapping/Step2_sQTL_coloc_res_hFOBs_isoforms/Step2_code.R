###########################################################################
###########################################################################
###                                                                     ###
###   sQTL coloc results and subset contained within hFOB isoforms      ###     
### 								                                                    ###
###########################################################################
###########################################################################

# get sQTLs that are in hFOBs 
##read in all sig events from all tissues and change them into a full dataframe
setwd("/Users/aa9gj/Documents/BPG_project/coloc_results_full/coloc_annot_results/")
my_files <- fread('my_files',header=FALSE)
my_data <- list()
for (i in seq_along(my_files$V1)) {
  my_data[[i]] <- fread(file = my_files$V1[i], header =TRUE)
  my_data[[i]]$tissue <- rep(my_files$V1[i], length(my_data[[i]]$V1))
}
coloc_sqtl_sig <- as.data.frame(do.call("rbind", my_data))
dim(coloc_sqtl_sig)#6889 total event identified by coloc sqtl
head(coloc_sqtl_sig)
coloc_sqtl_sig$event_id <- paste0(coloc_sqtl_sig$chr, "_", coloc_sqtl_sig$start, "_", coloc_sqtl_sig$end, "_", coloc_sqtl_sig$gene_name)
length(unique(coloc_sqtl_sig$event_id))#2043 unique events identified by sQTL coloc
sum(coloc_sqtl_sig$gene_type == "protein_coding")#6391
sum(coloc_sqtl_sig$gene_type != "protein_coding")#498
colnames(coloc_sqtl_sig) <- c('n_snps','H0','H1','H2','H3','H4','chr','start','end','clus','gene','gene_type','gene_name', 'tissue', 'event_id')
coloc_sqtl_sig$tissue[] <- sapply(coloc_sqtl_sig$tissue, gsub, pattern = "_annot", replacement = "")
coloc_sqtl_sig_pc <- filter(coloc_sqtl_sig, gene_type == 'protein_coding')#protein coding genes only
dim(coloc_sqtl_sig_pc)
str(coloc_sqtl_sig_pc)
length(unique(coloc_sqtl_sig_pc$event_id)) #1863 events with coloc sqtl in protein-coding genes
write.table(coloc_sqtl_sig_pc, file = "/Users/aa9gj/Documents/BPG_project/coloc_results_full/coloc_annot_results/coloc_sig_pc_events", quote = F, row.names = F)

#Group by gene, tissue, event
#for pc genes only
coloc_sqtl_sig_pc <- as.data.frame(coloc_sqtl_sig_pc)
by_gene_pc <- coloc_sqtl_sig_pc %>% group_by(gene_name, tissue) %>% dplyr::summarise(n = dplyr::n())
length(unique(by_gene_pc$gene_name))

gene_tissue_count <- table(by_gene_pc$gene_name)
gene_tissue_count_pc <- data.frame(gene_tissue_count)
write.table(gene_tissue_count_pc, "/Users/aa9gj/Documents/BPG_project/Materials_for_paper/Tables/S3", quote = F, row.names = F)

by_event_pc <- coloc_sqtl_sig_pc %>% group_by(event_id, tissue) %>% dplyr::summarise(n = n())
length(unique(by_event_pc$event_id))
by_tissue_pc <- coloc_sqtl_sig_pc %>% group_by(tissue) %>% dplyr::summarise(n = n())
length(unique(by_tissue_pc$tissue))
#write.table(by_tissue_pc, "/Users/aa9gj/Documents/BPG_project/Materials_for_paper/Tables/S2", quote = F, row.names = F)

# for all kinds of genes
by_gene <- coloc_sqtl_sig %>% group_by(gene_name, tissue) %>% dplyr::summarise(n = n())
length(unique(by_gene$gene_name))
by_gene_only <- coloc_sqtl_sig %>% group_by(gene_name) %>% dplyr::summarise(n = n())
length(unique(by_gene_only$gene_name))
by_event <- coloc_sqtl_sig %>% group_by(event_id, tissue) %>% dplyr::summarise(n = n())
length(unique(by_event$event_id))
by_tissue <- coloc_sqtl_sig %>% group_by(tissue) %>% dplyr::summarise(n = n())
length(unique(by_tissue$tissue))
#write.table(by_tissue, "/Users/aa9gj/Documents/BPG_project/Materials_for_paper/Tables/S1", quote = F, row.names = F)

by_event_summary <- coloc_sqtl_sig_pc %>% group_by(event_id) %>% dplyr::summarise(n = n())
by_event_summary <- do.call("data.frame", by_event_summary) # flatten
by_event_summary <- by_event_summary[order(by_event_summary$n, decreasing = T),]
by_event_summary$sep_col <- as.data.frame(str_split_fixed(by_event_summary$event_id, "_", 4))
by_event_summary <- do.call("data.frame",by_event_summary)
nrow(by_event_summary) #all protein coding genes with colocalizng sQTL
length(unique(by_event_summary$sep_col.V4))

## Perform analysis to figure out which of these sQTLs are within novel, known, or hybrid
events <- as.data.frame(by_event_summary$event_id)
colnames(events) <- "event_id"
events <- separate(events, "event_id", c("seq", "start", "end", "gene_name"), sep = "_")
hfob_coloc_genes <- as.data.frame(unique(by_event_summary$sep_col.V4))
colnames(hfob_coloc_genes) <- "gene_name"
## we have to add 1 to event start and take one off event end to match the introns constructed in sqanti output
events$start <- as.numeric(events$start) +1
events$end <- as.numeric(events$end) -1
events_gr <- makeGRangesFromDataFrame(events, keep.extra.columns = TRUE,
                                      start.field = "start", end.field = "end", seqnames.field = "seq", ignore.strand = T)

sqanti_gtf <- read_format("/Users/aa9gj/Documents/BPG_project/Materials_for_paper/aug2022_reanalysis/SQANTI3_results_full_corrected_chr_only.gtf")
sqanti_gtf_intron <- construct_introns(sqanti_gtf, update = TRUE)[]
sqanti_gtf_intron <- plyranges::filter(sqanti_gtf_intron, feature == "intron" | feature == "transcript")
sqanti_gtf_intron <- as.data.frame(sqanti_gtf_intron)


## Keep only isoforms belonging to these genes (whether they contain the event or not, and they are surviving the threshold)
isoforms_hfob_coloc <- inner_join(criteria1_3, hfob_coloc_genes, by = "gene_name")
length(unique(isoforms_hfob_coloc$gene_name))
isoforms_hfob_coloc <- inner_join(sqanti_gtf_intron, isoforms_hfob_coloc, by = c("transcript_id" = "rn"))
colnames(isoforms_hfob_coloc)
length(unique(isoforms_hfob_coloc$gene_name))
colnames(isoforms_hfob_coloc)
isoforms_hfob_coloc <- isoforms_hfob_coloc[,c(1:4, 6:10, 39:41, 90:91)]
isoforms_hfob_coloc_gr <- makeGRangesFromDataFrame(isoforms_hfob_coloc, keep.extra.columns = TRUE,
                                                   start.field = "start.x", end.field = "end.x", seqnames.field = "seqnames.x", ignore.strand = T, na.rm = TRUE)
length(unique(isoforms_hfob_coloc$gene_name))##coloc genes surviving the threshold without any exact match on events
length(unique(isoforms_hfob_coloc$transcript_id))##Number of transcripts belonging to those genes (all surviving threshold)
# Exact match overlap from start to end of event and intron
sqtl_overlap_isoforms_introns <- function(x,y) {
  hits  <- findOverlaps(query = x, subject = y, type = "equal")
  olap  <- pintersect(x[queryHits(hits)],
                      y[subjectHits(hits)])
  return(olap)
}

exact_overlap <- as.data.frame(sqtl_overlap_isoforms_introns(isoforms_hfob_coloc_gr, events_gr))#the isoform stays if it has the actual event from start to end
length(unique(exact_overlap$gene_name))#end up with 459
length(unique(exact_overlap$transcript_id))#2349

genes_in_hfobs_and_coloc <- as.data.frame(unique(exact_overlap$gene_name))
colnames(genes_in_hfobs_and_coloc) <- "gene_name"
length(unique(events$gene_name))
write.table(genes_in_hfobs_and_coloc, "/Users/aa9gj/Documents/BPG_project/Materials_for_paper/genes_in_hfobs_coloc", quote = F, row.names = F, col.names = F)
#make sure to add the +1 and -1 to get the events back to their original coordinates like in GTEx
exact_overlap$event_id <- paste0(exact_overlap$seqnames, "_", exact_overlap$start-1, "_", exact_overlap$end+1, "_", exact_overlap$gene_name)
length(unique(exact_overlap$event_id))#we end up with 836 unique events in coloc/hfobs/filter/exactmatch
# figure out the nature of those events
nature_of_events <- exact_overlap[, c("event_id", "transcript_id", "gene_name", "structural_category", "associated_transcript")]
#per event how many are in novel only, how many are in known only, and how many are in both
nature_of_events$status <- ifelse(nature_of_events$associated_transcript == "novel", "novel", "known")
summarize_events_s1 <- nature_of_events %>% group_by(event_id,status) %>% dplyr::summarize(n = dplyr::n())
summarize_events_s1 <- do.call("data.frame", summarize_events_s1) # flatten
summarize_events <- summarize_events_s1[,1:2]
summarize_events <- summarize_events %>% group_by(event_id) %>% dplyr::summarise(n = dplyr::n())
summarize_events <- do.call("data.frame", summarize_events) # flatten
filter(summarize_events, n == 2)#350 events explained by both known and novel isoforms 
test3 <- filter(summarize_events, n==1)
test3 <- inner_join(test3, summarize_events_s1, by = "event_id")
nrow(filter(test3, status == "known")) #383 events explained by known isoforms only
nrow(filter(test3, status == "novel")) #103 events explained by novel isoforms only

## Distribution of those genes per lead GWAS SNP, put it in a supplementary table
# create a genomic ranges for those genes 
dist_genes <- inner_join(genes_in_hfobs_and_coloc, gencodev38, by = "gene_name")
morris_lead_snps <- fread("/Users/aa9gj/Documents/BPG_project/morris_lead_snps_hg19.txt")
ch <- import.chain("/Users/aa9gj/Documents/BPG_project/hg19ToHg38.over.chain")
morris_overlift <- makeGRangesFromDataFrame(morris_lead_snps, keep.extra.columns = TRUE, start.field = "V4", end.field = "V4", seqnames.field = "V3", na.rm = T)
seqlevelsStyle(morris_overlift) = "UCSC"  # necessary
morris_overlift = liftOver(morris_overlift, ch)
class(morris_overlift)
morris_overlift = unlist(morris_overlift)
genome(morris_overlift) = "hg38"
morris_overlift <- as.data.frame(morris_overlift)
morris_overlift$snp_internal_id <- 1:nrow(morris_overlift)
morris_overlift <- morris_overlift[,-8]
morris_overlift_gr <- makeGRangesFromDataFrame(morris_overlift, keep.extra.columns = TRUE)
dist_genes <- dist_genes[,1:5]
dist_genes$start <- dist_genes$start - 200000
dist_genes$end <- dist_genes$end + 200000
dist_genes_gr <-  makeGRangesFromDataFrame(dist_genes, keep.extra.columns = TRUE)
# Exact match overlap from start to end of event and intron
dist_sgenes <- function(x,y) {
  hits  <- findOverlaps(query = x, subject = y, type = "any")
  olap  <- pintersect(x[queryHits(hits)],
                      y[subjectHits(hits)])
  return(olap)
}
test_dist1 <- as.data.frame(dist_sgenes(dist_genes_gr, morris_overlift_gr))
test_dist2 <- as.data.frame(dist_sgenes(morris_overlift_gr, dist_genes_gr))
write.table(test_dist1, "/Users/aa9gj/Documents/BPG_project/overlap_morris_sgenes1", quote = F, row.names = F)
write.table(test_dist2, "/Users/aa9gj/Documents/BPG_project/overlap_morris_sgenes2", quote = F, row.names = F)
length(unique(test_dist2$snp_internal_id))
length(unique(test_dist1$gene_name))
gwas_sgene_dist <- fread("/Users/aa9gj/Documents/BPG_project/gwas_sgenes_dist.txt")
gwas_sgene_dist_count <- gwas_sgene_dist %>% group_by(snp_internal_id) %>% dplyr::summarise(n = dplyr::n())
one_sgene <- filter(gwas_sgene_dist_count, n == 1)
more_sgene <- filter(gwas_sgene_dist_count, n > 1)
###########################################################################
###########################################################################
###                                                                     ###
###            stats on isoform containing the events                  ###     
### 								                                                    ###
###########################################################################
###########################################################################
## these isoforms have survived the threshold and the exact match, we can create another list that survived the threshold but didn't contain 
# the events for comparison purposes (later)
length(unique(exact_overlap$transcript_id))
length(unique(exact_overlap$gene_name))
novel_counts <- filter(exact_overlap, associated_transcript == "novel")
length(unique(novel_counts$transcript_id))
known_counts <- filter(exact_overlap, associated_transcript != "novel")
length(unique(known_counts$transcript_id))