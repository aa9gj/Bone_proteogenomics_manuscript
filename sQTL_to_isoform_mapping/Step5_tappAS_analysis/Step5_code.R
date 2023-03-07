###########################################################################
###########################################################################
###                                                                     ###
###                 DE, DIU, and their GO terms                         ###     
### 								                                                    ###
###########################################################################
###########################################################################
#set up gencode
gencode_v38 <- as.data.frame(rtracklayer::import('gencode.v38.annotation.gtf'))
gencode_v38_gene <- subset(gencode_v38, type == "gene")
gencode_v38_gene <- gencode_v38_gene[, c("gene_id", "gene_type", "gene_name")]
#keep only protein coding genes
gencode_v38_gene <- subset(gencode_v38_gene, gene_type == "protein_coding")
DEA_genes <- read.table("tappAS_DEA_Genes.tsv")
DEA_genes <- as.data.frame(left_join(DEA_genes, gencode_v38_gene, by = c("V1" = "gene_id")))
nrow(DEA_genes) #2034 DE genes in total
DE_genes <- filter(DEA_genes, DEA_genes$V3 == "DE")


DEA_transcript <- fread("tappAS_DEA_Transcripts.tsv")
DEA_transcript <- left_join(DEA_transcript, gencode_v38_gene, by = c("Gene" = "gene_id"))
gencode_v38_transcript <- subset(gencode_v38, type == "transcript")
gencode_v38_transcript <- gencode_v38_transcript[, 16:18]
gencode_v38_transcript <- subset(gencode_v38_transcript, transcript_type == "protein_coding")
DEA_transcript <- as.data.frame(left_join(DEA_transcript, gencode_v38_transcript, by = c("Name/Description" = "transcript_id")))
DE_transcript <- filter(DEA_transcript, DEA_transcript$`DEA Result` == "DE")


DIU <- fread("tappAS_DIUGene_Transcripts.tsv")
DIU <- as.data.frame(left_join(DIU, gencode_v38_gene, by = c("Gene Description" = "gene_id")))
DIU <- subset(DIU, DIU$`DIU Result` == "DIU") 
DIU_genes <- as.data.frame(DIU$gene_name)
colnames(DIU_genes) <- "gene_name"
