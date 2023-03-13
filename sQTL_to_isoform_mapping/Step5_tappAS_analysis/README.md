### All input files for tappAS analysis are found in the zenodo repo 

### Performing statistical analysis using [tappAS](https://tappas.org/)
1. Create a custom gff file from SQANTI output compatible with tappAS input using [IsoAnnotLite](https://isoannot.tappas.org/isoannot-lite/)
```shell
IsoAnnotLite_v2.7.0_SQ3.py SQANTI3_results_full_corrected.gtf SQANTI3_results_full_classification.txt SQANTI3_results_full_junctions.txt -gff3 GENCODE.gff3 -novel -o tappAS_full.gff3
```
2. Create a design file as input for tappAS depending on your needs. See tappAS webiste for further information. 
3. Create expression matrix file for input into tappAS (use cDNA_cupcake output)
4. Statistical analyses are performed using tappAS GUI as a time-series analysis using MaSigpro within tappAS (log fold change >= 2) for all DIU, DGE, and DTE
