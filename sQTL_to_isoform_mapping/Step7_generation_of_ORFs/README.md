### Where do I get these files to run this code? 
1. Download the following files from the Zenodo repo (coloc_annot_results; SQANTI3_results_full_corrected_chr_only.gtf, morris_lead_snps_hg19.txt)
2. Download GENCODE v38 gtf file from GENCODE
3. Download hg19ToHg38.over.chain from UCSC (chain file for liftOver)


### Calling open reading frames (ORFs)
We relied heavily on the already published [long-read proteogenomics (LRP) pipeline](https://github.com/sheynkman-lab/Long-Read-Proteogenomics) (with some modifications explained in the associated R code (Step7)

1. cpat_step (LRP)
```shell
cpat.py -x Human_Hexamer.tsv -d Human_logitModel.RData -g SQANTI3_results_full_corrected.fasta --min-orf=30 -o cpat_out_30minorf
```
2. orf_calling (LRP)
```shell
orf_calling_new.py --orf_coord cpat_out_30minorf.ORF_prob.tsv --orf_fasta cpat_out_30minorf.ORF_seqs.fa --gencode gencode.v38.annotation.gtf --sample_gtf SQANTI3_results_full_corrected.gtf --pb_gene pb_gene.tsv --classification SQANTI3_results_full_classification.txt --sample_fasta SQANTI3_results_full_corrected.fasta --output osteo_cpat.30minorf_ORF_called.tsv
```
3. refinement (LRP)
```shell
refine_orf_database.py --name hfobs --orfs osteo_cpat.30minorf_ORF_called.tsv --pb_fasta SQANTI3_results_full_corrected.fasta --coding_score_cutoff 0.0
```
4. make_cds (LRP)
```shell
make_cds_gtf_new.py --name hfob_CDS_with_transcripts_min30orf --sample_gtf SQANTI3_results_full_corrected.gtf --refined_database hfobs_orf_refined.tsv --called_orfs osteo_cpat.30minorf_ORF_called.tsv --pb_gene pb_gene.tsv --include_transcript yes
```
