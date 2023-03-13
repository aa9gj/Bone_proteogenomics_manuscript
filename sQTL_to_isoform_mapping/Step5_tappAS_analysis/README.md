### All input files for tappAS analysis are found in the zenodo repo 

### Performing statistical analysis using tappAS
1. Create a custom gff file from SQANTI output compatible with tappAS input using IsoAnnotLite (Shell script)
2. Create a design file as input for tappAS for each of the comparisons or time-series (R code)
3. Create expression matrix file for input into tappAS (R code)
4. Statistical analyses are performed using tappAS GUI as a time-series analysis using MaSigpro within tappAS (log fold change >= 2) for all DIU, DGE, and DTE
