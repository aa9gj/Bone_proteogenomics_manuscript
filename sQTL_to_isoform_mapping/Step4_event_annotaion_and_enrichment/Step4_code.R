###########################################################################
###########################################################################
###                                                                     ###
###                 events with effect sizes and their annotation       ###     
### 								                                                    ###
###########################################################################
###########################################################################

event_with_effect_size <- fread("/Users/aa9gj/Documents/BPG_project/Materials_for_paper/too_long_to_generate_files/event_with_effect_size")
event_with_effect_size <- separate(event_with_effect_size, phenotype_id, c("chr", "start", "end", "cluster", "gene_id"), sep = ":")
gencode_v26 <- as.data.frame(rtracklayer::import('/Users/aa9gj/Documents/BPG_project/Materials_for_paper/gencode.v26.annotation.gtf'))
gencode_v26 <- filter(gencode_v26, type == "gene", gene_type == "protein_coding")
gencode_v26 <- gencode_v26[,c("gene_id", "gene_name")]
event_with_effect_size <- left_join(event_with_effect_size, gencode_v26, by = "gene_id")
event_with_effect_size$event_id <- paste0(event_with_effect_size$chr, "_", event_with_effect_size$start, "_", event_with_effect_size$end, "_", event_with_effect_size$gene_name)
length(unique(event_with_effect_size$event_id))

events_hfob_coloc <- inner_join(exact_overlap, event_with_effect_size, by = "event_id")

# Now create the SNPs and their LD proxy (example here is for chr22, run a loop for all chromosomes)
chr22_lead_sqtl <- dplyr::filter(events_hfob_coloc, chr == "chr22")
chr22_bed <- fread("/Users/aa9gj/Documents/BPG_project/SNP_bed/chr22")
chr22_lead_sqtl_pos <- inner_join(chr22_bed, chr22_lead_sqtl, by = c("V4" = "V7"))
chr22_sqtl_and_LD <- list()
for (i in 1:nrow(chr22_lead_sqtl_pos)) {
  server <- "https://rest.ensembl.org"
  ext <- paste0("/ld/human/",chr22_lead_sqtl_pos$V4[i],"/1000GENOMES:phase_3:KHV?")
  r <- GET(paste(server, ext, sep = ""), content_type("application/json"))
  stop_for_status(r)
  chr22_sqtl_and_LD[[i]] <- fromJSON(toJSON(content(r)))
  print(i)
}
chr22_sqtl_and_LD_df <- as.data.frame(do.call("rbind", chr22_sqtl_and_LD))
myFun <- function(data) {
  ListCols <- sapply(data, is.list)
  cbind(data[!ListCols], t(apply(data[ListCols], 1, unlist)))
}
chr22_sqtl_and_LD_df <- myFun(chr22_sqtl_and_LD_df)
write.table(chr22_sqtl_and_LD_df, "/Users/aa9gj/Documents/BPG_project/chr_sqtl_ld/chr22_sqtl_and_ld_df", row.names = F, quote = F)
chr22_sqtl_and_LD_df <- filter(chr22_sqtl_and_LD_df, r2 >= 0.80)
chr22_sqtl_no_snps_LD <- as.data.frame(chr22_lead_sqtl_pos[which(!(chr22_lead_sqtl_pos$V4 %in% chr22_sqtl_and_LD_df$variation1)),])

## overlap with introns for 
chr22_lead_sqtl <- dplyr::filter(events_hfob_coloc, chr == "chr22")
chr22_sqtl_and_LD_df <- fread("/Users/aa9gj/Documents/BPG_project/chr_sqtl_ld/chr22_sqtl_and_ld_df")
chr22_sqtl_and_LD_df <- filter(chr22_sqtl_and_LD_df, r2 >= 0.80)
chr22_sqtl_no_snps_LD <- as.data.frame(chr22_lead_sqtl[which(!(chr22_lead_sqtl$V7 %in% chr22_sqtl_and_LD_df$variation1)),])
chr22_lead_no_LD <- as.data.frame(chr22_sqtl_no_snps_LD[,"V7"])
colnames(chr22_lead_no_LD) <- "rsid"
chr22_lead_no_LD <- as.data.frame(unique(chr22_lead_no_LD$rsid))
colnames(chr22_lead_no_LD) <- "rsid"

chr22_sqtl_LD_pos <- inner_join(chr22_bed, chr22_sqtl_and_LD_df, by = c("V4" = "variation2"))
chr22_lead_LD_pos <- chr22_sqtl_LD_pos[,"variation1"]
chr22_lead_LD_pos <- as.data.frame(unique(chr22_lead_LD_pos$variation1))
colnames(chr22_lead_LD_pos) <- "rsid"
chr22_lead_LD_pos <- right_join(chr22_bed, chr22_lead_LD_pos, by = c("V4" = "rsid"))
chr22_lead_LD_pos <- chr22_lead_LD_pos[,c(1:4)]
chr22_lead_LD_pos <- distinct(chr22_lead_LD_pos)

chr22_proxy_pos <- chr22_sqtl_LD_pos[,c("V1", "V2", "V3", "V4", "r2", "d_prime", "variation1")]
chr22_proxy_pos <- distinct(chr22_proxy_pos)

chr22_lead_no_LD_pos <- right_join(chr22_bed, chr22_lead_no_LD, by = c("V4" = "rsid"))
chr22_lead_no_LD_pos <- chr22_lead_no_LD_pos[,c(1:4)]
chr22_lead_no_LD_pos <- distinct(chr22_lead_no_LD_pos)

event_location <- as.data.frame(unique(chr22_lead_sqtl$event_id))
colnames(event_location) <- "event"
event_location <- separate(event_location, "event", c("seq", "start", "end", "gene"), sep = "_")
event_location$seq <- gsub("chr", "", event_location$seq)
event_location_gr <- makeGRangesFromDataFrame(event_location, keep.extra.columns = TRUE,
                                              start.field = "start", end.field = "end", seqnames.field = "seq", ignore.strand = T)
chr22_lead_no_LD_pos_gr <- makeGRangesFromDataFrame(chr22_lead_no_LD_pos, keep.extra.columns = TRUE,
                                                    start.field = "V2", end.field = "V3", seqnames.field = "V1", ignore.strand = T, na.rm = TRUE)

chr22_lead_LD_pos_gr <- makeGRangesFromDataFrame(chr22_lead_LD_pos, keep.extra.columns = TRUE,
                                                 start.field = "V2", end.field = "V3", seqnames.field = "V1", ignore.strand = T)

sqtl_overlap_introns <- function(x,y) {
  hits  <- findOverlaps(query = x, subject = y, type = "any")
  olap  <- pintersect(x[queryHits(hits)],
                      y[subjectHits(hits)])
  return(olap)
}

chr22_noLD_lead_overlap <- as.data.frame(sqtl_overlap_introns(chr22_lead_no_LD_pos_gr, event_location_gr))
chr22_LD_lead_overlap <- as.data.frame(sqtl_overlap_introns(chr22_lead_LD_pos_gr, event_location_gr))

chr22_proxy_wo_overlap_lead <- subset(chr22_proxy_pos, !(chr22_proxy_pos$variation1 %in% chr22_LD_lead_overlap$V4))
chr22_proxy_wo_overlap_lead_gr <- makeGRangesFromDataFrame(chr22_proxy_wo_overlap_lead, keep.extra.columns = TRUE,
                                                           start.field = "V2", end.field = "V3", seqnames.field = "V1", ignore.strand = T)
chr22_proxy_overlap <- as.data.frame(sqtl_overlap_introns(chr22_proxy_wo_overlap_lead_gr, event_location_gr))


##create introns from sqanti gtf file for ssa, ssd overlaps

sqanti_gtf <- read_format("/Users/aa9gj/Documents/BPG_project/Materials_for_paper/aug2022_reanalysis//SQANTI3_results_full_corrected_chr_only.gtf")
isoforms_wanted <- as.data.frame(exact_overlap[, "transcript_id"])
colnames(isoforms_wanted) <- "transcript_id"
sqanti_gtf_intron <- construct_introns(sqanti_gtf, update = TRUE)[]
sqanti_gtf_intron <- plyranges::filter(sqanti_gtf_intron, feature == "intron")
sqanti_gtf_intron_threshold <- inner_join(isoforms_wanted, sqanti_gtf_intron, by = "transcript_id", copy = T)

sqanti_gtf_intron_threshold$ssd_start <- sqanti_gtf_intron_threshold$start
sqanti_gtf_intron_threshold$ssd_end <- sqanti_gtf_intron_threshold$start+1
sqanti_gtf_intron_threshold$ssa_start <- sqanti_gtf_intron_threshold$end-2
sqanti_gtf_intron_threshold$ssa_end <- sqanti_gtf_intron_threshold$end-1
sqanti_gtf_intron_threshold$seqnames <- gsub("chr", "", sqanti_gtf_intron_threshold$seqnames)
ssd_coord <- sqanti_gtf_intron_threshold[,c(1,2,6,12,13)]
ssa_coord <- sqanti_gtf_intron_threshold[,c(1,2,6,14,15)]
ssd_coord_gr <- makeGRangesFromDataFrame(ssd_coord, keep.extra.columns = TRUE,
                                         start.field = "ssd_start", end.field = "ssd_end", seqnames.field = "seqnames", ignore.strand = F)
ssa_coord_gr <- makeGRangesFromDataFrame(ssa_coord, keep.extra.columns = TRUE,
                                         start.field = "ssa_start", end.field = "ssa_end", seqnames.field = "seqnames", ignore.strand = F)

chr22_lead_sqtl <- dplyr::filter(events_hfob_coloc, chr == "chr22")
chr22_sqtl_and_LD_df <- fread("/Users/aa9gj/Documents/BPG_project/chr_sqtl_ld/chr22_sqtl_and_ld_df")
chr22_sqtl_and_LD_df <- filter(chr22_sqtl_and_LD_df, r2 >= 0.80)
chr22_sqtl_no_snps_LD <- as.data.frame(chr22_lead_sqtl[which(!(chr22_lead_sqtl$V7 %in% chr22_sqtl_and_LD_df$variation1)),])
chr22_lead_no_LD <- as.data.frame(chr22_sqtl_no_snps_LD[,"V7"])
colnames(chr22_lead_no_LD) <- "rsid"
chr22_lead_no_LD <- as.data.frame(unique(chr22_lead_no_LD$rsid))
colnames(chr22_lead_no_LD) <- "rsid"
chr22_bed <- fread("/Users/aa9gj/Documents/BPG_project/SNP_bed/chr22")

chr22_sqtl_LD_pos <- inner_join(chr22_bed, chr22_sqtl_and_LD_df, by = c("V4" = "variation2"))
chr22_lead_LD_pos <- chr22_sqtl_LD_pos[,"variation1"]
chr22_lead_LD_pos <- as.data.frame(unique(chr22_lead_LD_pos$variation1))
colnames(chr22_lead_LD_pos) <- "rsid"
chr22_lead_LD_pos <- right_join(chr22_bed, chr22_lead_LD_pos, by = c("V4" = "rsid"))
chr22_lead_LD_pos <- chr22_lead_LD_pos[,c(1:4)]
chr22_lead_LD_pos <- distinct(chr22_lead_LD_pos)

chr22_proxy_pos <- chr22_sqtl_LD_pos[,c("V1", "V2", "V3", "V4", "r2", "d_prime", "variation1")]
chr22_proxy_pos <- distinct(chr22_proxy_pos)

chr22_lead_no_LD_pos <- right_join(chr22_bed, chr22_lead_no_LD, by = c("V4" = "rsid"))
chr22_lead_no_LD_pos <- chr22_lead_no_LD_pos[,c(1:4)]
chr22_lead_no_LD_pos <- distinct(chr22_lead_no_LD_pos)

chr22_lead_no_LD_pos_gr <- makeGRangesFromDataFrame(chr22_lead_no_LD_pos, keep.extra.columns = TRUE,
                                                    start.field = "V2", end.field = "V3", seqnames.field = "V1", ignore.strand = T, na.rm = TRUE)

chr22_lead_LD_pos_gr <- makeGRangesFromDataFrame(chr22_lead_LD_pos, keep.extra.columns = TRUE,
                                                 start.field = "V2", end.field = "V3", seqnames.field = "V1", ignore.strand = T)


chr22_proxy_pos_gr <- makeGRangesFromDataFrame(chr22_proxy_pos, keep.extra.columns = TRUE,
                                               start.field = "V2", end.field = "V3", seqnames.field = "V1", ignore.strand = T)

sqtl_overlap_splice_sites <- function(x,y) {
  hits  <- findOverlaps(query = x, subject = y, type = "any")
  olap  <- pintersect(x[queryHits(hits)],
                      y[subjectHits(hits)])
  return(olap)
}

ssd_hits_part1 <- sqtl_overlap_splice_sites(chr22_lead_no_LD_pos_gr, ssd_coord_gr)
ssd_hits_part2 <- sqtl_overlap_splice_sites(chr22_lead_LD_pos_gr, ssd_coord_gr)
ssd_hits_part3 <- sqtl_overlap_splice_sites(chr22_proxy_pos_gr, ssd_coord_gr)

ssa_hits_part1 <- sqtl_overlap_splice_sites(chr22_lead_no_LD_pos_gr, ssa_coord_gr)
ssa_hits_part2 <- sqtl_overlap_splice_sites(chr22_lead_LD_pos_gr, ssa_coord_gr)
ssa_hits_part3 <- sqtl_overlap_splice_sites(chr22_proxy_pos_gr, ssa_coord_gr)

ssd_hits_part1
ssd_hits_part2
ssd_hits_part3

ssa_hits_part1
ssa_hits_part2
ssa_hits_part3

### # let's do enrichment analysis 
## get the lead sQTL and their LD and their positions
chr22_lead_sqtl <- dplyr::filter(events_hfob_coloc, chr == "chr22")
chr22_sqtl_and_LD_df <- fread("/Users/aa9gj/Documents/BPG_project/chr_sqtl_ld/chr22_sqtl_and_ld_df")
chr22_sqtl_and_LD_df <- filter(chr22_sqtl_and_LD_df, r2 >= 0.80)
chr22_sqtl_no_snps_LD <- as.data.frame(chr22_lead_sqtl[which(!(chr22_lead_sqtl$V7 %in% chr22_sqtl_and_LD_df$variation1)),])
chr22_lead_no_LD <- as.data.frame(chr22_sqtl_no_snps_LD[,"V7"])
colnames(chr22_lead_no_LD) <- "rsid"
chr22_lead_no_LD <- as.data.frame(unique(chr22_lead_no_LD$rsid))
colnames(chr22_lead_no_LD) <- "rsid"

chr22_bed <- fread("/Users/aa9gj/Documents/BPG_project/SNP_bed/chr22")
chr22_sqtl_LD_pos <- inner_join(chr22_bed, chr22_sqtl_and_LD_df, by = c("V4" = "variation2"))
chr22_lead_LD_pos <- chr22_sqtl_LD_pos[,"variation1"]
chr22_lead_LD_pos <- as.data.frame(unique(chr22_lead_LD_pos$variation1))
colnames(chr22_lead_LD_pos) <- "rsid"
chr22_lead_LD_pos <- right_join(chr22_bed, chr22_lead_LD_pos, by = c("V4" = "rsid"))
chr22_lead_LD_pos <- chr22_lead_LD_pos[,c(1:4)]
chr22_lead_LD_pos <- distinct(chr22_lead_LD_pos)

chr22_proxy_pos <- chr22_sqtl_LD_pos[,c("V1", "V2", "V3", "V4", "r2", "d_prime", "variation1")]
chr22_proxy_pos <- distinct(chr22_proxy_pos)

chr22_lead_no_LD_pos <- right_join(chr22_bed, chr22_lead_no_LD, by = c("V4" = "rsid"))
chr22_lead_no_LD_pos <- chr22_lead_no_LD_pos[,c(1:4)]
chr22_lead_no_LD_pos <- distinct(chr22_lead_no_LD_pos)

chr22 <- rbind(chr22_lead_LD_pos, chr22_lead_no_LD_pos)
chr22 <- rbind(chr22, chr22_proxy_pos[,1:4])
rm(chr22_bed)
rm(chr22_lead_LD_pos)
rm(chr22_lead_no_LD_pos)
write.table(chr22, "/Users/aa9gj/Documents/BPG_project/test_enrichment_in_r/lead_ld_position_by_chr/chr22", row.names = F, quote = F)

setwd("/Users/aa9gj/Documents/BPG_project/test_enrichment_in_r/lead_ld_position_by_chr/")
my_files <- fread('chromosomes',header=FALSE)
my_data <- list()
for (i in seq_along(my_files$V1)) {
  my_data[[i]] <- fread(file = my_files$V1[i], header =TRUE)
}
lead_sqtl_and_proxy <- as.data.frame(do.call("rbind", my_data))

lead_sqtl_and_proxy <- distinct(lead_sqtl_and_proxy)

sqtl_overlap_SF_binding <- function(x,y) {
  hits  <- findOverlaps(query = x, subject = y, type = "any")
  olap  <- pintersect(x[queryHits(hits)],
                      y[subjectHits(hits)])
  return(olap)
}

library(vsea)
vsea_lead_sqtl_full <- as.data.frame(unique(events_hfob_coloc$V7))
colnames(vsea_lead_sqtl_full) <- "rsid"
lead_sqtl_annotated <- inner_join(vsea_lead_sqtl_full, lead_sqtl_and_proxy, by = c('rsid' = "V4"))
length(unique(lead_sqtl_annotated$rsid))
ch <- import.chain("/Users/aa9gj/Documents/BPG_project/test_enrichment_in_r/hg38ToHg19.over.chain")
cur <- makeGRangesFromDataFrame(lead_sqtl_annotated, keep.extra.columns = TRUE, start.field = "V3", end.field = "V3", seqnames.field = "V1", na.rm = T)
seqlevelsStyle(cur) = "UCSC"  # necessary
cur19 = liftOver(cur, ch)
class(cur19)
cur19 = unlist(cur19)
genome(cur19) = "hg19"
cur19_df <- as.data.frame(cur19)
cur19_df$seqnames <- gsub("chr", "", cur19_df$seqnames)
cur19_df <- cur19_df[,c(1,3)]
cur19_df <- paste0(cur19_df$seqnames, ":", cur19_df$end)
snpsnap <- read_snpsnap(population = c("EUR"))
bg_variants <- select_background_variants(cur19_df, snpsnap, bg_size = 1)
bg_variants <- bg_variant[!is.na(bg_variants)]
bg_variants <- as.data.frame(do.call(rbind, bg_variants))
bg_variants <- as.data.frame(bg_variants[,1])
colnames(bg_variants) <- "snp"
bg_variants$chr <- rep("chr", nrow(bg_variants), )
bg_variants <- separate(bg_variants, "snp", c("seqname", "start"))
bg_variants$end <- bg_variants$start
bg_variants$seqnames <- paste0(bg_variants$chr, bg_variants$seqname)
bg_variants_gr <- makeGRangesFromDataFrame(bg_variants, keep.extra.columns = TRUE, start.field = "start", end.field = "end", seqnames.field = "seqnames", na.rm = T)
#automate_sf_bind
meta_eclip <- fread("/Users/aa9gj/Documents/BPG_project/test_enrichment_in_r/metadata.tsv")
meta_eclip <- meta_eclip[, c(1,3,11, 23)]
meta_eclip$`Experiment target` <- gsub("-human", "", meta_eclip$`Experiment target`)
colnames(meta_eclip) <- c("eclip_id", "type", "celltype" ,"gene")
## how many SF in our samples
yeo_nature <- fread("/Users/aa9gj/Documents/BPG_project/Actual_SFs_from_yeo.txt", header = F)
hfob_coloc_yeo_all <- genes_in_hfobs_and_coloc[which(genes_in_hfobs_and_coloc$gene_name %in% yeo_nature$V1),]
yeo_eclip <- inner_join(meta_eclip, yeo_nature, by = c("gene" = "V1"))
length(unique(yeo_eclip$gene))
yeo_eclip <- filter(yeo_eclip, type == "bed")
length(unique(yeo_eclip$gene))
eclip_sf_ids <- as.data.frame(yeo_eclip$eclip_id)
colnames(eclip_sf_ids) <- "id"
list_enrich <- list()
for (i in 1:nrow(eclip_sf_ids)) {
  splice_factor <- fread(paste0("/Users/aa9gj/Documents/BPG_project/test_enrichment_in_r/liftedhg38tohg19_bed_eclip/",eclip_sf_ids[i,],"_lift.bed"))
  splice_factor_gr <- makeGRangesFromDataFrame(splice_factor, keep.extra.columns = TRUE,
                                               start.field = "V2", end.field = "V3", seqnames.field = "V1", ignore.strand = T)
  sqtl_set <- sqtl_overlap_SF_binding(cur19, splice_factor_gr)
  test_set <- sqtl_overlap_SF_binding(test_no_na_gr, splice_factor_gr)
  lead_value <- nrow(as.data.frame(sqtl_set))
  back_value <- nrow(as.data.frame(test_set))
  
  dat <- data.frame(
    "within" = c(lead_value, back_value),
    "out" = c(1572-lead_value, 1572-back_value),
    row.names = c("lead", "back"),
    stringsAsFactors = FALSE
  )
  colnames(dat) <- c("within", "out")
  
  SF_enrich_results <- fisher.test(dat)
  list_enrich[[i]] <- paste(SF_enrich_results$p.value,",", yeo_eclip$gene[i], "," ,yeo_eclip$celltype[i])
}
enrich_df <- as.data.frame(do.call("rbind", list_enrich))
enrich_df <- separate(enrich_df, "V1", c("p_value", "SF", "cell"), sep = ",")
sig_sf <- filter(enrich_df, p_value < 0.05)
length(unique(sig_sf$SF))
length(unique(enrich_df$SF))

enrich_df$fdr <- p.adjust(enrich_df$p_value, method = "fdr", n = nrow(enrich_df))
write.table(enrich_df, "/Users/aa9gj/Documents/BPG_project/test_enrichment_in_r/enrich_df", row.names = F, quote = F)
enrich_df <- fread("/Users/aa9gj/Documents/BPG_project/test_enrichment_in_r/enrich_df")
sig_sf_FDR <- filter(enrich_df, fdr < 0.05)
length(unique(sig_sf_FDR$SF))
in_coloc_hfob <- as.data.frame(unique(events_hfob_coloc$gene_name.y))
colnames(in_coloc_hfob) <- "gene"
enrich_df$SF <- gsub(" ", "", enrich_df$SF)
in_coloc_hfob_enrich <- inner_join(enrich_df, in_coloc_hfob, by = c("SF" = "gene"))