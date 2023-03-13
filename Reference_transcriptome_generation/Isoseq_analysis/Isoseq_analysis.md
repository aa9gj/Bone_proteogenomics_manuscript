### Performing long-read RNAseq analysis from raw reads to isoforms classification 

1. Raw reads can be found in GEO (ref number) please allocate enough space for download. 
2. Primer removal and demultiplexing analysis using [lima](https://lima.how/) and the list of barcodes provided in (file name)
3. Refine using isoseq refine module to trim of poly(A) tails and perform rapid concatemer identification and removal
4. Cluster using isoseq cluster module using hierarchical n*log(n) alignment and iterative cluster merging and then Polished POA sequence generation, using a QV guided consensus approach
5. Align reads to the reference transcriptome using isoseq align
6. Perform isoform classification using [SQANTI3](https://github.com/ConesaLab/SQANTI3) 
