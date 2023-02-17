# generation of source data for prioritization scheme
#set up gencode
gencode_v38 <- as.data.frame(rtracklayer::import('/Users/aa9gj/Documents/BPG_project/Materials_for_paper/gencode.v38.annotation.gtf'))
gencode_v38_gene <- subset(gencode_v38, type == "gene")
gencode_v38_gene <- gencode_v38_gene[, c("gene_id", "gene_type", "gene_name")]
#keep only protein coding genes
gencode_v38_gene <- subset(gencode_v38_gene, gene_type == "protein_coding")
DEA_genes <- read.table("/Users/aa9gj/Documents/BPG_project/Materials_for_paper/tappAS_oct8th2022/tappAS_DEA_Genes.tsv")
DEA_genes <- as.data.frame(left_join(DEA_genes, gencode_v38_gene, by = c("V1" = "gene_id")))
DE_genes <- filter(DEA_genes, DEA_genes$V3 == "DE")


DEA_transcript <- fread("/Users/aa9gj/Documents/BPG_project/Materials_for_paper/tappAS_oct8th2022/tappAS_DEA_Transcripts.tsv")
DEA_transcript <- left_join(DEA_transcript, gencode_v38_gene, by = c("Gene" = "gene_id"))
gencode_v38_transcript <- subset(gencode_v38, type == "transcript")
gencode_v38_transcript <- gencode_v38_transcript[, 16:18]
gencode_v38_transcript <- subset(gencode_v38_transcript, transcript_type == "protein_coding")
DEA_transcript <- as.data.frame(left_join(DEA_transcript, gencode_v38_transcript, by = c("Name/Description" = "transcript_id")))
DE_transcript <- filter(DEA_transcript, DEA_transcript$`DEA Result` == "DE")

DIU <- fread("/Users/aa9gj/Documents/BPG_project/Materials_for_paper/tappAS_oct8th2022/tappAS_DIUGene_Transcripts.tsv")
DIU <- as.data.frame(left_join(DIU, gencode_v38_gene, by = c("Gene Description" = "gene_id")))
DIU <- subset(DIU, DIU$`DIU Result` == "DIU") 

## any in IMPC?
abnormal_bmd <- fread("/Users/aa9gj/Documents/BPG_project/Materials_for_paper/supporting_data/abnormal_bone_mineral_density.tsv", header = T)
abnormal_bmd_genes <- as.data.frame(toupper(abnormal_bmd$Gene))
colnames(abnormal_bmd_genes) <- "gene_name"
increased_bmd <- fread("/Users/aa9gj/Documents/BPG_project/Materials_for_paper/supporting_data/increased_bone_mineral_density.tsv", header = T)
increased_bmd_genes <- as.data.frame(toupper(increased_bmd$Gene))
colnames(increased_bmd_genes) <- "gene_name"
decreased_bmd <- fread("/Users/aa9gj/Documents/BPG_project/Materials_for_paper/supporting_data/decreased_bone_mineral_density.tsv", header = T)
decreased_bmd_genes <- as.data.frame(toupper(decreased_bmd$Gene))
colnames(decreased_bmd_genes) <- "gene_name"
bmd_genes <- full_join(abnormal_bmd_genes, increased_bmd_genes, by = "gene_name")
bmd_genes <- full_join(bmd_genes, decreased_bmd_genes, by = "gene_name")
bmd_genes <- as.data.frame(unique(bmd_genes$gene_name))
colnames(bmd_genes) <- "gene_name"
bmd_genes$gene_name <- gsub("ZFP", "ZNF", bmd_genes$gene_name)

## any in bone literature?
bans_genes <- fread("/Users/aa9gj/Documents/BPG_project/Materials_for_paper/supporting_data/Basel_BANs.csv", header = T)
bans_genes <- bans_genes[BAN == "YES",]
bans_genes <- as.data.frame(toupper(bans_genes$Gene))
colnames(bans_genes) <- "gene_name"
bans_genes$gene_name <- gsub("ZFP", "ZNF", bans_genes$gene_name)

eqtl_genes <- fread("/Users/aa9gj/Documents/BPG_project/Materials_for_paper/supporting_data/Basel_eqtl.txt", header = T)
eqtl_genes <- as.data.frame(unique(eqtl_genes$`Gene Name`))
colnames(eqtl_genes) <- "gene_name"

mono_genes <- fread("/Users/aa9gj/Documents/BPG_project/Materials_for_paper/supporting_data/bone_genes_olivia", header = F)
mono_genes <- as.data.frame(toupper(mono_genes$V1))
colnames(mono_genes) <- "gene_name"
mono_genes$gene_name <- gsub("ZFP", "ZNF", mono_genes$gene_name)

# Create source data for the 459 genes 
source_data_tableS <- as.data.frame(unique(exact_overlap$event_id))
colnames(source_data_tableS) <- "event_id"
source_data_tableS <- separate(source_data_tableS, "event_id", c("chr", "start", "end", "gene_name"))
source_data_tableS$event_id <- paste0(source_data_tableS$chr, "_", source_data_tableS$start, "_", source_data_tableS$end, "_", source_data_tableS$gene_name)
source_data_tableS$DE_gene <- ifelse(source_data_tableS$gene_name %in% DEA_genes$gene_name, 1, 0)
exact_overlap_source <- exact_overlap[, c(11,12,14,18)]
exact_overlap_source <- distinct(exact_overlap_source)
source_data_tableS <- inner_join(exact_overlap_source, source_data_tableS, by = "event_id")
source_data_tableS$DE_transcript <- ifelse(source_data_tableS$transcript_id %in% DE_transcript$`#Transcript`, 1, 0)
source_data_tableS$DIU <- ifelse(source_data_tableS$gene_name %in% DIU$gene_name, 1, 0)  
major_iso_switch <- filter(DIU, DIU$`Major Isoform Switching` == "YES")
source_data_tableS$major_isoform_switch <- ifelse(source_data_tableS$gene_name %in% major_iso_switch$gene_name, 1, 0)  
source_data_tableS$impc_evidence <- ifelse(source_data_tableS$gene_name %in% bmd_genes$gene_name, 1, 0)  
source_data_tableS$BAN_evidence <- ifelse(source_data_tableS$gene_name %in% bans_genes$gene_name, 1, 0)
source_data_tableS$eQTL_evidence <- ifelse(source_data_tableS$gene_name %in% eqtl_genes$gene_name, 1, 0)
source_data_tableS$monogenic_evidence <- ifelse(source_data_tableS$gene_name %in% mono_genes$gene_name, 1, 0)
source_data_tableS <- inner_join(source_data_tableS, gene_tissue_count_pc, by = c("gene_name" = "Var1"))
source_data_tableS <- inner_join(source_data_tableS, by_event_pc_source, by = c("event_id" = "Var1"))
write.table(source_data_tableS, "/Users/aa9gj/Documents/BPG_project/Materials_for_paper/source_data_feb6", quote = F, row.names = F)
