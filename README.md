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
  - **sQTL_colocalization_analysis**: This directory contains code needed to replicate Bayesian colocalization analysis with [Coloc](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1004383). Please refer to the README.md within directory for further information
    - Step 0: Perform bayesian colocalization analysis using summary statistics from the latest [BMD GWAS](https://www.nature.com/articles/s41588-018-0302-x) with summary statistics from sQTL data for all 49 GTEx tissues. 
  - **Reference_transcriptome_generation**: This directory contains code to generate the reference transcriptome from long-read RNAseq data. Please refer to the README.md within directory for further information
    - Step 1: Perform analyses on outputs from SQANTI and cDNA_cupcake
  - **sQTL_to_isoform_mapping**
    - Step 2: Characterize full-length isoforms (known and novel) containing the colocalized junctions
    - Step 3: description
    - Step 4: description
    - Step 5: description
    - Step 6: description
    - Step 7: description
    - Step 8: NMD and truncation analysis was performed using a beta version of [Biosurfer](https://github.com/sheynkman-lab/Biosurfer_BMD_analysis) 
