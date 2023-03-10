## Code accompanying the manuscript *"Long read proteogenomics to connect disease-associated sQTLs to the protein isoform effectors in disease"* 

The full text can be found in [Abood et al. 2023, BioRxiv](https://www.biorxiv.org/) 


The data required is found in DOI 10.5281/zenodo.7603851 <br>

The raw long-read sequencing data is found in [GSE224588](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE224588)

### How to use this repository

1. Use setup_r_env.R to set up the R environment with all the needed packages. 
2. The repo is broken down into three major sections: 
  - **sQTL_colocalization_analysis**: This repo contains code needed to replicate Bayesian colocalization analysis with [Coloc](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1004383)
  - **Reference_transcriptome_generation**
  - **sQTL_to_isoform_mapping**
