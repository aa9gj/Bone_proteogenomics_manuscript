echo "track name=\"BMD coloc sQTL junctions\" description=\"BMD coloc sQTL junctions - from Arby\"" > tmp.bed
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $5, $6, $4, $7}' coloc_res.bed > coloc_res_correct_columns.bed
cat tmp.bed coloc_res_correct_columns.bed > coloc_res_correct_columns_final.bed
rm tmp.bed
