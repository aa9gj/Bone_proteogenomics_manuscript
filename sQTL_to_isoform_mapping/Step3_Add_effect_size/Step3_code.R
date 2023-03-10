# Where do I get these files to run this code? 
# Download the following files from the Zenodo repo (Morris_data/chr22.morris.gwas.snps; Gtex_lookup/chr22.gtex.lookup.txt, event_snps/chr22.events)

# Adding effect size from both GTEx summary stats and GWAS summary stats. This should be done in a loop and takes a long time.  
# Ensure that you run this in a loop, this is just an example for chromosome 22
morris <- fread("Morris_data/chr22.morris.gwas.snps")
prep_gwas <- function(x) {
  x$V1 <- gsub("chr", "", x$V1)
  x<-x[order(x$V16, decreasing = F),]
  return(x)
}
gwas <- prep_gwas(morris)
colnames(gwas) <- c("seqnames", "start", "end", "width", "strand", "SNPID", "rsid", "EA","NEA","EAF","INFO","BETA","SE","P.nominal","P.I","P.bmd","N","maf","ID")
gtex_lookup <- fread("Gtex_lookup/chr22.gtex.lookup.txt")

# Get the lowest and then the second lowest rsid and pvalue and find out the difference
prep_sqtl <- function(x,y) {
  colnames(x[[i]]) <- c("phenotype_id", "variant_id", "tss_dist", "ma_samples", "ma_count", "maf", "pval_nominal", "slope", "slope_se")
  x[[i]] <- inner_join(x[[i]], y, by = c("variant_id"="V1"))
  x[[i]]<-x[[i]][order(x[[i]]$pval_nominal, decreasing = F),]
  return(x[[i]])
}

my_snps <- list()
event_files <- fread("all_events_snps/chr22.events", header = F)
setwd("all_events_snps")
for (i in seq_along(event_files$V1)) {
  my_snps[[i]] <- fread(event_files$V1[i])
  my_snps[[i]] <- prep_sqtl(my_snps, gtex_lookup)
  my_snps[[i]]$effect_size <- gwas[grep(my_snps[[i]]$V7[1], gwas$rsid),]$BETA[1]
  my_snps[[i]]$pval_diff <- my_snps[[i]]$pval_nominal[1] - my_snps[[i]]$pval_nominal[2]
  my_snps[[i]] <- my_snps[[i]][1, c(1,6:9,12,13,15,17,18)]
  print(i)
}

chr22_events <- as.data.frame(do.call("rbind", my_snps))
# This data is fully available on Zenodo so you don't have to run the code above to get it but it is there for your consideration. 
write.table(chr22_events, "chr22_events_effect", row.names = F, quote = F)

# Once you generate all chr??_events you can run this
eventslist <- list()
event_files <- fread("chr_events_for_effect_sizes/event_files", header = F)
setwd("chr_events_for_effect_sizes/sQTL/")
for (i in seq_along(event_files$V1)) {
  eventslist[[i]] <- fread(event_files$V1[i])
}
