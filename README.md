[![DOI](https://zenodo.org/badge/597128359.svg)](https://zenodo.org/badge/latestdoi/597128359)
## Code accompanying the manuscript *"Long read proteogenomics to connect disease-associated sQTLs to the protein isoform effectors in disease"* 

The full text can be found in [Abood et al. 2023, BioRxiv](https://www.biorxiv.org/content/10.1101/2023.03.17.531557v1) 

### Purpose 

We present a novel generalizable approach that integrates information from GWAS, splicing QTL (sQTL), and PacBio long-read RNA-seq in a disease relevant model to infer the effects of sQTLs on the ultimate protein isoform products they encode

### Data availability 
1. Processed and input data is found in [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7603851.svg)](https://doi.org/10.5281/zenodo.7603851)
2. Raw long-read sequencing data is found in [GSE224588](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE224588)

### How to use this repository

1. Use setup_r_env.R to set up the R environment with all the needed packages. 
2. The repo is broken down into three major sections: 
  - **sQTL_colocalization_analysis**: This directory contains code needed to replicate Bayesian colocalization analysis with [Coloc](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1004383). Please refer to the README.md within directory for further information
    - Step 0: Perform bayesian colocalization analysis using summary statistics from the latest [BMD GWAS](https://www.nature.com/articles/s41588-018-0302-x) with summary statistics from sQTL data for all 49 GTEx tissues. 
  - **Reference_transcriptome_generation**: This directory contains code to generate the reference transcriptome from long-read RNAseq data. Please refer to the README.md within directory for further information
    - [Isoseq analysis](https://github.com/aa9gj/Bone_proteogenomics_manuscript/blob/main/Reference_transcriptome_generation/Isoseq_analysis/Isoseq_analysis.md): from raw reads to isoform classification
    - [Step 1](https://github.com/aa9gj/Bone_proteogenomics_manuscript/tree/main/Reference_transcriptome_generation/Step1_Long-read_RNAseq_filtering_in_hFOBs): Perform analyses on outputs from SQANTI and cDNA_cupcake
  - **sQTL_to_isoform_mapping**
    - Step 2: Characterize full-length isoforms (known and novel) containing the colocalized junctions
    - Step 3: Add effect size and direction of effect to colocalized junctions
    - Step 4: Annotate lead sQTLs and their proxy, follow with positional and enrichment analyses
    - Step 5: Differential analyses (DE and DIU) using tappAS
    - Step 6: Integrating multiple datasets from the literature and within our analyses to prioritize the isoforms for experimental validation
    - Step 7: ORF analyses including: NMD and truncation analysis was performed using a beta version of [Biosurfer](https://github.com/sheynkman-lab/Biosurfer_BMD_analysis) 
