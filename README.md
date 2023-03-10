## Code accompanying the manuscript *"Long read proteogenomics to connect disease-associated sQTLs to the protein isoform effectors in disease"* 

The full text can be found in [Abood et al. 2023, BioRxiv](https://www.biorxiv.org/) 

### Purpose 

We present a novel generalizable approach that integrates information from GWAS, splicing QTL (sQTL), and PacBio long-read RNA-seq in a disease relevant model to infer the effects of sQTLs on the ultimate protein isoform products they encode

### Data availability 
Processed data is found in DOI 10.5281/zenodo.7603851 <br>
Raw long-read sequencing data is found in [GSE224588](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE224588)

### How to use this repository

1. Use setup_r_env.R to set up the R environment with all the needed packages. 
2. The repo is broken down into three major sections: 
  - **sQTL_colocalization_analysis**: This directory contains code needed to replicate Bayesian colocalization analysis with [Coloc](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1004383)
    - Step 0: Please refer to the README.md within directory for further information
  - **Reference_transcriptome_generation**
    - Step 1: Please refer to the README.md within directory for further information
  - **sQTL_to_isoform_mapping**
    - Step 2: description
    - Step 3: description
    - Step 4: 
    - Step 5: 
    - Step 6:
    - Step 7: 
    - Step 8: NMD and truncation analysis was performed using a Beta version of [Biosurfer](https://github.com/sheynkman-lab/Biosurfer_BMD_analysis) 
