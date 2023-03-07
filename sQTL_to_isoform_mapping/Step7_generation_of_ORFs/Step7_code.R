###########################################################################
###########################################################################
###                                                                     ###
###                LRP predicted CDS hfob (needs to be updated as of March 7th                               ###     
### 								                                                    ###
###########################################################################
###########################################################################
LRP_gtf <- as.data.frame(rtracklayer::import("/Users/aa9gj/Documents/BPG_project/hfob_CDS_with_transcripts_min30orf_with_cds.gtf"))


# Do this without getting rid of the exons
LRP_gtf <- separate(LRP_gtf, "transcript_id", c("gene_name", "transcript_id", "extra_column"), sep = "\\|")
LRP_gtf <- LRP_gtf[, c(1:7, 11:12)]
colnames(LRP_gtf) <- c("seqnames", "CDS_start", "CDS_end", "CDS_width", "CDS_strand", "CDS_source","CDS_type", "CDS_gene_name", "transcript_id")
LRP_gtf_annot <- inner_join(criteria1_3, LRP_gtf, by = c("gene_name" = "CDS_gene_name"))
length(unique(LRP_gtf_annot$gene_name))
LRP_gtf_annot_coloc <- inner_join(exact_overlap, LRP_gtf, by = c("gene_name" = "CDS_gene_name"))
length(unique(LRP_gtf_annot_coloc$gene_name))

refined_cds <- fread("/Users/aa9gj/Documents/BPG_project/Biosurfer_input/hfobs_att2_orf_refined.tsv")
fix_test1 <- inner_join(refined_cds, LRP_gtf_annot_coloc, by = c("base_acc" = "transcript_id.y"))
length(unique(fix_test1$transcript_id))
#filter(fix_test1, CDS_type == "transcript")
target_transcripts <- as.data.frame(unique(exact_overlap$transcript_id))
colnames(target_transcripts) <- "transcript_id"

new_gtf <- list()
for (i in seq_along(target_transcripts$transcript_id)) {
  new_gtf[[i]] <- fix_test1[grep(target_transcripts$transcript_id[i], fix_test1$pb_accs),]
  new_gtf[[i]]$target_transcript <- target_transcripts$transcript_id[i]
}

new_gtf_df <- do.call("rbind", new_gtf)
length(unique(new_gtf_df$target_transcript))#2045 transcripts
length(unique(new_gtf_df$gene_name))#445 genes

new_gtf_df <- new_gtf_df[,c(21:25, 27:35)]
colnames(new_gtf_df)[7] <- "seqnames"
colnames(new_gtf_df)[14] <- "transcript_id"

export(new_gtf_df, "/Users/aa9gj/Desktop/CDS.gtf", format = 'gtf')
test_gtf <- as.data.frame(rtracklayer::import("/Users/aa9gj/Desktop/CDS.gtf"))

# overlap events on whehter they are within coding sequence or UTR 
test_gtf <- test_gtf[,-c(4,5)]
test_gtf_gr <- makeGRangesFromDataFrame(test_gtf, keep.extra.columns = TRUE,
                                        start.field = "start", end.field = "end", seqnames.field = "seqnames", ignore.strand = T)

overlap_all <- sqtl_overlap_UTR_CDS(test_gtf_gr,events_gr)
overlap_all_df <- as.data.frame(overlap_all)
length(unique(overlap_all_df$transcript_id))
length(unique(overlap_all_df$gene_name))

test_gtf_CDS_only <- filter(test_gtf, CDS_type == "CDS")
test_gtf_CDS_only <- test_gtf_CDS_only[,-c(4,5)]
test_gtf_UTR_only <- filter(test_gtf, CDS_type == "exon")
test_gtf_UTR_only <- test_gtf_UTR_only[,-c(4,5)]

test_gtf_CDS_only <- makeGRangesFromDataFrame(test_gtf_CDS_only, keep.extra.columns = TRUE,
                                              start.field = "start", end.field = "end", seqnames.field = "seqnames", ignore.strand = T)
test_gtf_UTR_only <- makeGRangesFromDataFrame(test_gtf_UTR_only, keep.extra.columns = TRUE,
                                              start.field = "start", end.field = "end", seqnames.field = "seqnames", ignore.strand = T)

sqtl_overlap_UTR_CDS <- function(x,y) {
  hits  <- findOverlaps(query = x, subject = y, type = "any")
  olap  <- pintersect(x[queryHits(hits)],
                      y[subjectHits(hits)])
  return(olap)
}

length(unique(test_gtf$transcript_id))
length(unique(test_gtf$gene_name))
overlap_CDS <- sqtl_overlap_UTR_CDS(test_gtf_CDS_only,events_gr)
overlap_CDS_df <- as.data.frame(overlap_CDS)
length(unique(overlap_CDS_df$transcript_id))
length(unique(overlap_CDS_df$gene_name))

overlap_CDS_df[grep("PTBP1", overlap_CDS_df$gene_name),]#Its there

# Create the gtf file and support information for Mayank 
trans_abund_mayank <- as.data.frame(unique(test_gtf$transcript_id))
colnames(trans_abund_mayank) <- "transcript_id"
trans_abund_mayank <- inner_join(trans_abund_mayank, criteria1_3, by = c("transcript_id" = "rn"))
trans_abund_mayank <- trans_abund_mayank[,1:12]
write.table(trans_abund_mayank, "/Users/aa9gj/Documents/BPG_project/Materials_for_paper/aug2022_reanalysis/Biosurfer_input_for_Mayank/trans_abund", quote = F, row.names = F)

step1 <- as.data.frame(unique(test_gtf$gene_name))
colnames(step1) <- "gene_name"
step1 <- inner_join(step1, criteria1_3, by = "gene_name")
trans_abund_myank_not <- step1[which(!(step1$rn %in% trans_abund_mayank$transcript_id)),]
trans_abund_myank_not <- trans_abund_myank_not[,2:13]
write.table(trans_abund_myank_not, "/Users/aa9gj/Documents/BPG_project/Materials_for_paper/aug2022_reanalysis/Biosurfer_input_for_Mayank/trans_abund_not", quote = F, row.names = F)

# create the non-overlap gtf file for mayank and now he has everything 
step1 <- as.data.frame(unique(test_gtf$gene_name))
colnames(step1) <- "gene_name"
step1 <- inner_join(step1, criteria1_3, by = "gene_name")
trans_abund_myank_not <- step1[which(!(step1$rn %in% trans_abund_mayank$transcript_id)),]
trans_abund_myank_not <- trans_abund_myank_not[,1:13]
not_exact_gene <- as.data.frame(unique(trans_abund_myank_not$gene_name))
colnames(not_exact_gene) <- "gene_name"
LRP_gtf_annot_coloc <- inner_join(not_exact_gene_trans, LRP_gtf, by = c("gene_name" = "CDS_gene_name"))
length(unique(LRP_gtf_annot_coloc$gene_name))

refined_cds <- fread("/Users/aa9gj/Documents/BPG_project/Biosurfer_input/hfobs_att2_orf_refined.tsv")
fix_test1 <- inner_join(refined_cds, LRP_gtf_annot_coloc, by = c("base_acc" = "transcript_id"))
length(unique(fix_test1$pb_accs))
#filter(fix_test1, CDS_type == "transcript")
target_transcripts <- as.data.frame(unique(trans_abund_myank_not$rn))
colnames(target_transcripts) <- "transcript_id"

new_gtf_not_exact <- list()
for (i in seq_along(target_transcripts$transcript_id)) {
  new_gtf_not_exact[[i]] <- fix_test1[grep(target_transcripts$transcript_id[i], fix_test1$pb_accs),]
  new_gtf_not_exact[[i]]$target_transcript <- target_transcripts$transcript_id[i]
}

new_gtf_not_exact_df <- do.call("rbind", new_gtf_not_exact)
length(unique(new_gtf_not_exact_df$target_transcript))#2045 transcripts
length(unique(new_gtf_not_exact_df$gene_name))#445 genes

new_gtf_not_exact_df <- new_gtf_not_exact_df[,c(10:18)]
colnames(new_gtf_df)[7] <- "seqnames"
colnames(new_gtf_df)[14] <- "transcript_id"


export(new_gtf_not_exact_df, "/Users/aa9gj/Desktop/CDS_not_exact.gtf", format = 'gtf')
test_gtf <- as.data.frame(rtracklayer::import("/Users/aa9gj/Desktop/CDS.gtf"))
