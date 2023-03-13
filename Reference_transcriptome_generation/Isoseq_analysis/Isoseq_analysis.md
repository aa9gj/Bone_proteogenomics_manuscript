### Performing long-read RNAseq analysis from raw reads to isoforms classification 

1. Raw reads can be found in GEO ([GSE224588](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE224588)). Please allocate enough space for download. 
2. Primer removal and demultiplexing analysis using [lima](https://lima.how/) and the list of barcodes provided in (file name)
3. Refine using isoseq refine module to trim of poly(A) tails and perform rapid concatemer identification and removal
4. Cluster using isoseq cluster module using hierarchical n*log(n) alignment and iterative cluster merging and then Polished POA sequence generation, using a QV guided consensus approach
5. Generate raw isoforms counts post clustering using cDNA_cupcake module (demux_isoseq_with_genome.py)
```shell
ml anaconda/2020.11-py3.8
export PATH=$PATH:/cDNA_Cupcake/sequence/
export PATH=$PATH:/cDNA_Cupcake/rarefaction/

demux_isoseq_with_genome.py --mapped_fafq --read_stat --classify_csv merged_flnc_report -o counts_cupcake.mapped_fl_count.txt
```
/home/aa9gj/bin/cDNA_Cupcake/post_isoseq_cluster/demux_isoseq_with_genome.py --mapped_fafq new_all_transcriptome.fasta --read_stat new_all_transcriptome.read_stat.txt --classify_csv merged_flnc_report -o counts_cupcake.mapped_fl_count.txt
7. Align reads to the reference transcriptome using isoseq align
8. Perform isoform classification using [SQANTI3](https://github.com/ConesaLab/SQANTI3) 
