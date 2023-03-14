###########################################################################
###########################################################################
###                                                                     ###
###                LRP predicted CDS hfob                               ###     
### 								                                                    ###
###########################################################################
###########################################################################
# Read the LRP output (GTF file)
LRP_gtf <- as.data.frame(rtracklayer::import("hfob_CDS_with_transcripts_min30orf_with_cds.gtf"))
LRP_gtf <- separate(LRP_gtf, "transcript_id", c("gene_name", "transcript_id", "extra_column"), sep = "\\|")
LRP_gtf <- LRP_gtf[, c(1:7, 11:12)]
colnames(LRP_gtf) <- c("seqnames", "CDS_start", "CDS_end", "CDS_width", "CDS_strand", "CDS_source","CDS_type", "CDS_gene_name", "CDS_transcript_id")
# Keep ORFs that have passed the long-read expression threshold (supplemental note 2) 
LRP_gtf_annot <- inner_join(criteria1_3, LRP_gtf, by = c("gene_name" = "CDS_gene_name"))
# Keep ORFs that have an exact match sQTL
LRP_gtf_annot_coloc <- inner_join(exact_overlap, LRP_gtf, by = c("gene_name" = "CDS_gene_name"))
# Read in the refined CDS file from the LRP pipeline
refined_cds <- fread("hfobs_att2_orf_refined.tsv")
test_genes <- as.data.frame(unique(criteria1_3$gene_name))

LRP_fixed <- inner_join(refined_cds, LRP_gtf_annot_coloc, by = c("base_acc" = "CDS_transcript_id"))
target_transcripts <- as.data.frame(unique(exact_overlap$transcript_id))
colnames(target_transcripts) <- "transcript_id"

new_gtf <- list()
for (i in seq_along(target_transcripts$transcript_id)) {
  new_gtf[[i]] <- LRP_fixed[grep(paste0("\\b",target_transcripts$transcript_id[i],"\\b"), fix_test1$pb_accs),]
  new_gtf[[i]]$target_transcript <- target_transcripts$transcript_id[i]
}

# Here is the corrected LRP output gtf file
LRP_fixed <- do.call("rbind", new_gtf)


