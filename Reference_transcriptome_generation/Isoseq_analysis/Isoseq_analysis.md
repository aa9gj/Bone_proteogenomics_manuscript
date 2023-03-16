### Performing long-read RNAseq analysis from raw reads to isoforms classification 

1. Raw reads can be found in GEO ([GSE224588](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE224588)). Please allocate enough space for download. 

2. Refine using isoseq [refine module](https://isoseq.how/getting-started.html) to trim of poly(A) tails and perform rapid concatemer identification and removal
```shell
isoseq3 refine $lima_ouput_file barcodes.fasta output_file.flnc.bam --require-polya -j 8
```
3. Cluster using isoseq [cluster module](https://isoseq.how/getting-started.html) using hierarchical n*log(n) alignment and iterative cluster merging and then Polished POA sequence generation, using a QV guided consensus approach
```shell
isoseq3 cluster flnc.fofn clustered.bam --verbose --use-qvs -j 40
```
3. Generate raw isoforms counts post clustering using [cDNA_cupcake](https://github.com/Magdoll/cDNA_Cupcake) module demux_isoseq_with_genome.py
```shell
ml anaconda/2020.11-py3.8
export PATH=$PATH:/cDNA_Cupcake/sequence/
export PATH=$PATH:/cDNA_Cupcake/rarefaction/

demux_isoseq_with_genome.py --mapped_fafq --read_stat --classify_csv merged_flnc_report -o counts_cupcake.mapped_fl_count.txt
```
5. Align reads to the reference transcriptome using [minimap2](https://github.com/PacificBiosciences/pbmm2)
```shell
pbmm2 align clustered.bam GRCh38.primary_assembly.genome.mmi aligned_sorted.bam --preset ISOSEQ --sort
```
6. Collapse redundant transcripts (based on exonic structures) into unique isoforms using isoseq3 [collapse module](https://isoseq.how/classification/isoseq-collapse.html)
```shell
isoseq3 collapse aligned_sorted.bam new_all_transcriptome.gff
```
7. Perform isoform classification using [SQANTI3](https://github.com/ConesaLab/SQANTI3) 
```shell
sqanti3_qc.py --fasta .fasta gencode.v38.annotation.gtf GRCh38.primary_assembly.genome.fa --skipORF -o SQANTI3_results_full --fl transcriptome.abundance.txt
```
