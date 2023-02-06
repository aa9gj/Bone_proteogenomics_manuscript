#!/bin/bash

#SBATCH -J make_gene_files
#SBATCH -n 1
#SBATCH -t 72:00:00
#SBATCH -o out
#SBATCH -e err
#SBATCH -p standard
#SBATCH -A cphg-farber
#SBATCH --mem=100000

filename='list_of_directories'
filelines=`cat ${filename}`
#echo Start
for line in ${filelines} ; do
    (cd $line && ls *.genes > ./$line.file_names && cd ..)
done
	
