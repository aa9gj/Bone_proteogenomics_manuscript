# adding effect size 

morris <- fread("/Users/aa9gj/Documents/BPG_project/Mirrorplots_data_code_results/Morris_data/morris_snps_bychr/chr22.morris.gwas.snps")
prep_gwas <- function(x) {
  x$V1 <- gsub("chr", "", x$V1)
  x<-x[order(x$V16, decreasing = F),]
  return(x)
}
gwas <- prep_gwas(morris)
colnames(gwas) <- c("seqnames", "start", "end", "width", "strand", "SNPID", "rsid", "EA","NEA","EAF","INFO","BETA","SE","P.nominal","P.I","P.bmd","N","maf","ID")


gtex_lookup <- fread("/Users/aa9gj/Documents/BPG_project/Mirrorplots_data_code_results/Gtex_lookup/GTEx_lookup_table_bychr/chr22.gtex.lookup.txt")

# get the lowest and then the second lowest rsid and pvalue
prep_sqtl <- function(x,y) {
  colnames(x[[i]]) <- c("phenotype_id", "variant_id", "tss_dist", "ma_samples", "ma_count", "maf", "pval_nominal", "slope", "slope_se")
  x[[i]] <- inner_join(x[[i]], y, by = c("variant_id"="V1"))
  x[[i]]<-x[[i]][order(x[[i]]$pval_nominal, decreasing = F),]
  return(x[[i]])
}

my_snps <- list()
event_files <- fread("/Users/aa9gj/Documents/BPG_project/Mirrorplots_data_code_results/event_snps/event_with_sig_colocalization/chr22.events", header = F)
setwd("/Users/aa9gj/Documents/BPG_project/Mirrorplots_data_code_results/event_snps/event_with_sig_colocalization/all_events_snps")
for (i in seq_along(event_files$V1)) {
  my_snps[[i]] <- fread(event_files$V1[i])
  my_snps[[i]] <- prep_sqtl(my_snps, gtex_lookup)
  my_snps[[i]]$effect_size <- gwas[grep(my_snps[[i]]$V7[1], gwas$rsid),]$BETA[1]
  my_snps[[i]]$pval_diff <- my_snps[[i]]$pval_nominal[1] - my_snps[[i]]$pval_nominal[2]
  my_snps[[i]] <- my_snps[[i]][1, c(1,6:9,12,13,15,17,18)]
  print(i)
}

#write this as a loop to get all the chromosomes (see example chr1)
chr1_events <- as.data.frame(do.call("rbind", my_snps))

write.table(chr1_events, "/Users/aa9gj/Documents/BPG_project/chr_events_for_effect_sizes/sQTL/chr1_events_effect", row.names = F, quote = F)

eventslist <- list()
event_files <- fread("/Users/aa9gj/Documents/BPG_project/chr_events_for_effect_sizes/event_files", header = F)
setwd("/Users/aa9gj/Documents/BPG_project/chr_events_for_effect_sizes/sQTL/")
for (i in seq_along(event_files$V1)) {
  eventslist[[i]] <- fread(event_files$V1[i])
}