#!/bin/bash

i=$SLURM_ARRAY_TASK_ID
mine=$(sed -n "${i}{p;}" < sqtl_output_dirs_names.txt)
fine=$(sed -n "${i}{p;}" < full_file_names.txt)

filename='chr1_'
filelines=`cat ${filename}`
#echo Start
for line in ${filelines} ; do
    zgrep $line ${fine}  > ${mine}/$line.genes
done
