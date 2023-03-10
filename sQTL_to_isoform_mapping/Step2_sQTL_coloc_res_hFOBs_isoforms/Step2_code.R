###########################################################################
###########################################################################
###                                                                     ###
###   sQTL coloc results and subset contained within hFOB isoforms      ###     
### 								                                                    ###
###########################################################################
###########################################################################

# Where do I get these files to run this code? 
# Download the following files from the Zenodo repo (coloc_annot_results; SQANTI3_results_full_corrected_chr_only.gtf, morris_lead_snps_hg19.txt)
# Download GENCODE v38 gtf file from GENCODE
# Download hg19ToHg38.over.chain from UCSC (chain file for liftOver)

# Get sQTLs that are in hFOBs 
# Read in all junctions with colocalized sQTL (significant) from all tissues into a list and convert them into a dataframe
setwd("coloc_annot_results/")
my_files <- fread('my_files',header=FALSE)
my_data <- list()
for (i in seq_along(my_files$V1)) {
  my_data[[i]] <- fread(file = my_files$V1[i], header =TRUE)
  my_data[[i]]$tissue <- rep(my_files$V1[i], length(my_data[[i]]$V1))
}

# Junctions with colocalizing sQTLs in all gene types
coloc_sqtl_sig <- as.data.frame(do.call("rbind", my_data))
dim(coloc_sqtl_sig)#6889 total event identified by coloc sqtl
coloc_sqtl_sig$event_id <- paste0(coloc_sqtl_sig$chr, "_", coloc_sqtl_sig$start, "_", coloc_sqtl_sig$end, "_", coloc_sqtl_sig$gene_name)
length(unique(coloc_sqtl_sig$event_id))#2043 unique events identified by sQTL coloc

# Create a dataframe with significant junctions with colocalizing sQTLs 
colnames(coloc_sqtl_sig) <- c('n_snps','H0','H1','H2','H3','H4','chr','start','end','clus','gene','gene_type','gene_name', 'tissue', 'event_id')
coloc_sqtl_sig$tissue[] <- sapply(coloc_sqtl_sig$tissue, gsub, pattern = "_annot", replacement = "")
coloc_sqtl_sig_pc <- filter(coloc_sqtl_sig, gene_type == 'protein_coding')#protein coding genes only
length(unique(coloc_sqtl_sig_pc$event_id)) #1863 events with coloc sqtl in protein-coding genes

by_event_summary <- coloc_sqtl_sig_pc %>% group_by(event_id) %>% dplyr::summarise(n = n())
by_event_summary <- do.call("data.frame", by_event_summary) # flatten
by_event_summary <- by_event_summary[order(by_event_summary$n, decreasing = T),]
by_event_summary$sep_col <- as.data.frame(str_split_fixed(by_event_summary$event_id, "_", 4))
by_event_summary <- do.call("data.frame",by_event_summary)
nrow(by_event_summary) #all protein coding genes with colocalizng sQTL
length(unique(by_event_summary$sep_col.V4))

# Perform analysis to figure out which of these sQTLs are within novel, known, or hybrid
events <- as.data.frame(by_event_summary$event_id)
colnames(events) <- "event_id"
events <- separate(events, "event_id", c("seq", "start", "end", "gene_name"), sep = "_")
hfob_coloc_genes <- as.data.frame(unique(by_event_summary$sep_col.V4))
colnames(hfob_coloc_genes) <- "gene_name"

# Add 1 to event start and take one off event end to match the introns constructed in sqanti output
events$start <- as.numeric(events$start) +1
events$end <- as.numeric(events$end) -1
events_gr <- makeGRangesFromDataFrame(events, keep.extra.columns = TRUE,
                                      start.field = "start", end.field = "end", seqnames.field = "seq", ignore.strand = T)

sqanti_gtf <- read_format("SQANTI3_results_full_corrected_chr_only.gtf")
sqanti_gtf_intron <- construct_introns(sqanti_gtf, update = TRUE)[]
sqanti_gtf_intron <- plyranges::filter(sqanti_gtf_intron, feature == "intron" | feature == "transcript")
sqanti_gtf_intron <- as.data.frame(sqanti_gtf_intron)

# Keep only isoforms belonging to these genes (whether they contain the event or not, and they are surviving the threshold)
isoforms_hfob_coloc <- inner_join(criteria1_3, hfob_coloc_genes, by = "gene_name")
isoforms_hfob_coloc <- inner_join(sqanti_gtf_intron, isoforms_hfob_coloc, by = c("transcript_id" = "rn"))
isoforms_hfob_coloc <- isoforms_hfob_coloc[,c(1:4, 6:10, 39:41, 90:91)]
isoforms_hfob_coloc_gr <- makeGRangesFromDataFrame(isoforms_hfob_coloc, keep.extra.columns = TRUE,
                                                   start.field = "start.x", end.field = "end.x", seqnames.field = "seqnames.x", ignore.strand = T, na.rm = TRUE)

# Exact match overlap from start to end of event and intron
sqtl_overlap_isoforms_introns <- function(x,y) {
  hits  <- findOverlaps(query = x, subject = y, type = "equal")
  olap  <- pintersect(x[queryHits(hits)],
                      y[subjectHits(hits)])
  return(olap)
}

# The isoform stays if it has the actual event from start to end
exact_overlap <- as.data.frame(sqtl_overlap_isoforms_introns(isoforms_hfob_coloc_gr, events_gr))
length(unique(exact_overlap$gene_name))#end up with 459 genes
length(unique(exact_overlap$transcript_id))# with 2349 isoforms
genes_in_hfobs_and_coloc <- as.data.frame(unique(exact_overlap$gene_name))
colnames(genes_in_hfobs_and_coloc) <- "gene_name"

# Make sure to add the +1 and -1 to get the events back to their original coordinates like in GTEx
exact_overlap$event_id <- paste0(exact_overlap$seqnames, "_", exact_overlap$start-1, "_", exact_overlap$end+1, "_", exact_overlap$gene_name)
length(unique(exact_overlap$event_id))#we end up with 836 unique events in coloc/hfobs/filter/exactmatch

# Figure out the nature of those events
nature_of_events <- exact_overlap[, c("event_id", "transcript_id", "gene_name", "structural_category", "associated_transcript")]
# Per event how many are in novel only, how many are in known only, and how many are in both
nature_of_events$status <- ifelse(nature_of_events$associated_transcript == "novel", "novel", "known")
summarize_events_s1 <- nature_of_events %>% group_by(event_id,status) %>% dplyr::summarize(n = dplyr::n())
summarize_events_s1 <- do.call("data.frame", summarize_events_s1) # flatten
summarize_events <- summarize_events_s1[,1:2]
summarize_events <- summarize_events %>% group_by(event_id) %>% dplyr::summarise(n = dplyr::n())
summarize_events <- do.call("data.frame", summarize_events) # flatten
filter(summarize_events, n == 2) #350 events explained by both known and novel isoforms 
nat_events <- filter(summarize_events, n==1)
nat_events <- inner_join(nat_events, summarize_events_s1, by = "event_id")
nrow(filter(nat_events, status == "known")) #383 events explained by known isoforms only
nrow(filter(nat_events, status == "novel")) #103 events explained by novel isoforms only within hFOB

# Distribution of those genes per lead GWAS SNP, put it in a supplementary table
# Create a genomic ranges for those genes and keep things within 400 Kb
dist_genes <- inner_join(genes_in_hfobs_and_coloc, gencodev38, by = "gene_name")
morris_lead_snps <- fread("morris_lead_snps_hg19.txt")
ch <- import.chain("hg19ToHg38.over.chain")
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
write.table(test_dist1, "overlap_morris_sgenes1", quote = F, row.names = F)
write.table(test_dist2, "overlap_morris_sgenes2", quote = F, row.names = F)
######## upload gwas_sgenes_dist.txt to zenodo
gwas_sgene_dist <- fread("gwas_sgenes_dist.txt")
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
# These isoforms have survived the threshold and the exact match

novel <- filter(exact_overlap, associated_transcript == "novel")
known <- filter(exact_overlap, associated_transcript != "novel")
